import os
import subprocess
import ringdown
from . import File
import pandas as pd
import numpy as np
import h5py

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from .File import *
from .StrainDatabase import *
from .PosteriorDatabase import *
from . import File
from . import StrainDatabase
from . import PosteriorDatabase

from . import metadb

class Database:
    def __init__(self, data_folder, posterior_urls=None, strain_urls=None, psd_urls=None):
        if data_folder is not None:
            if data_folder[-1] == "/":
                data_folder
            self.data_folder = data_folder
            self.posterior_folder = f"{data_folder}/PosteriorData"
            self.strain_folder = f"{data_folder}/StrainData"
        else:
            self.data_folder = None
            self.posterior_folder = None
            self.strain_folder = None

        # Pull the default url data for all the posteriors, strains and psds. 
        self.posterior_urls = None
        if posterior_urls is None:
            posterior_url_path = pkg_resources.open_text(metadb, 'posterior_urls.csv')
            self.posterior_urls = pd.read_csv(posterior_url_path)

        self.strain_urls = None
        if strain_urls is None:
            strain_url_path = pkg_resources.open_text(metadb, 'strain_urls.csv')
            self.strain_urls = pd.read_csv(strain_url_path)

        self.psd_urls = None
        if psd_urls is None:
            psd_url_path = pkg_resources.open_text(metadb, 'psd_urls.csv')
            self.psd_urls = pd.read_csv(psd_url_path)

    def initialize(self, data_folder=None):
        # This will overwrite the default folder if folder is provided:
        if data_folder is not None:
            if data_folder[-1] == "/":
                data_folder = data_folder[0:-1]
            self.data_folder = data_folder
            self.posterior_folder = f"{data_folder}/PosteriorData"
            self.strain_folder = f"{data_folder}/StrainData"

        # Create the folder if it doesn't exist, and then create the
        # posterior and strain folders as well
        for folder in [self.data_folder, self.posterior_folder, self.strain_folder]:
            if not os.path.exists(folder):
                print(f"Creating a folder at {folder}")
                subprocess.run(["mkdir", folder])

        # This will create the databases
        self.PosteriorDB = PosteriorDatabase(self.posterior_folder, self.posterior_urls, self.psd_urls, self.strain_urls)
        self.StrainDB = StrainDatabase(self.strain_folder, self.strain_urls)

    def update_posterior_schema(self, schema_addition):
        self.PosteriorDB.schema.update(schema_addition)

    def update_strain_schema(self, schema_addition):
        self.StrainDB.schema.update(schema_addition)

    def event(self, eventname):
        return Event(eventname, self)

    def event_list(self):
        return list(self.strain_urls.Event.unique())


class Event:
    def __init__(self, eventname, DB_reference):
        self.name = eventname
        self.DB_ref = DB_reference
        self.PD_ref = DB_reference.PosteriorDB
        self.SD_ref = DB_reference.StrainDB

    def posteriors(self):
        """
        Returns a dataframe of the posterior samples of the event

        Returns:
            pd.DataFrame: Posterior samples of the event with each row
            being a posterior sample and the columns the parameter
        """
        return self.PD_ref.posteriors(self.name)

    def psd(self, detector=None):
        """
        Returns the PSD (Power Spectral Density) estimate for all detectors
        or a single detector if specified.
        
        Will download the needed PSD files if they are not already there

        Args:
            detector (None, string, or list of strings):
                If not provided, it will return the psds of all detectors 
                available.
                A list of detectors whose PSDs you want can be input as
                e.g. detector = ['H1','L1']
                A string denoting the detector whose psd you want 'H1',
                'L1' or 'V1'

        Returns:
            A dictionary containing ringdown.PowerSpectrum objects for
            each detector's PSD.

            If detector is not specified:
                A dictionary labelled by detector name, containing PSDs
                for all detectors as ringdown.PowerSpectrum objects

            If detector is specified as a string (e.g. detector='H1'):
                A ringdown.PowerSpectrum object containing the detector PSD
        """
        return self.PD_ref.psd(self.name, detector=detector)

    def strain(self, detectors=None, duration=32.0):
        """
        Returns the strain for all detectors or a single detector 
        if specified.
        
        Will download the needed hdf5 files if they are not already there

        Args:
            detector (None, string, or list of strings):
                If not provided, it will return the strain object of all detectors 
                available.
                A list of detectors whose PSDs you want can be input as
                e.g. detector = ['H1','L1']
                A string denoting the detector whose psd you want 'H1',
                'L1' or 'V1'

            duration (float):
                default is 32.0s. You can ask for a 4096.0s strain as well
                but as of now it probably will not overwrite what you have.
                This is definitely on the next TODO.

        Returns:
            A dictionary containing ringdown.PowerSpectrum objects for
            each detector's PSD.

            If detector is not specified:
                A dictionary labelled by detector name, containing PSDs
                for all detectors as ringdown.PowerSpectrum objects

            If detector is specified as a string (e.g. detector='H1'):
                A ringdown.PowerSpectrum object containing the detector PSD
        """
        return self.SD_ref.strain(self.name, detectors=detectors, duration=duration)

    def read_posterior_file(self, h5path, datatype='array', attr_name=None, detectors=None, approximant=None, replacement_dict=None):
        """
        Returns the referenced data object in the posterior hdf5 file for all detectors or a single detector 
        if specified.

        Will not work for GWTC-1 events, and will only work once the file 
        is already downloaded

        Args:
            h5path (string):
                Internal path for data you want from a particular hdf5 file.
                Use '{detector}', '{approximant}', '{event}' as wildcards to 
                be filled in with the appropriate values for this event. 
                If {detector} is used as a wildcard, results for each appropriate
                detector will be returned.

            datatype (string):
                Can be one of ('array', 'value', 'attribute').
                If it is an attribute, you must also provide attr_name='name-of-attribute'

        Returns:
            A dictionary containing the asked for objects for
            each detector's PSD.

            If detector is not specified:
                A dictionary labelled by detector name, containing objects 
                for each detector

            If detector is specified as a string (e.g. detector='H1'):
                returns the object for that one detector value
        """
        # Run over the list of detectors if not provided
        if ((detectors is None) and ("{detector}" in h5path)):
            detectors = self.PD_ref.available_detectors(self.name)

        # Choose the best approximant according to the provided list order
        if approximant is None:
            approximant = self.PD_ref.choose_approximant(self.name)

        # Method to interpolate the path
        if replacement_dict is None:
            replacement_dict = {'detector': detectors,  'event': self.name, 'approximant': approximant}

        # if the detectors are a list, return a list of dictionaries
        if isinstance(detectors, list):
            result = {}
            for ifo in detectors:
                scheme = {'type': datatype, 'name': attr_name, 'path': h5path}
                file = self.PD_ref.event_path(self.name)
                new_replacement_dict = replacement_dict.copy()
                new_replacement_dict['detector'] = ifo
                with h5py.File(file, 'r') as f:
                    result[ifo] = self.PD_ref.read_data_from_file(f, scheme, new_replacement_dict)
        else:
            scheme = {'type': datatype, 'name': attr_name, 'path': h5path}
            file = self.PD_ref.event_path(self.name)
            with h5py.File(file, 'r') as f:
                print(scheme)
                print(replacement_dict)
                result = self.PD_ref.read_data_from_file(f, scheme, replacement_dict)
        return result

    def read_strain_file(self, h5path, datatype='array', attr_name=None, detectors=None, replacement_dict=None):
        """
        Returns the referenced data object in the strain hdf5 file for all detectors or a single detector 
        if specified.

        Will not work for GWTC-1 events, and will only work once the file 
        is already downloaded

        Args:
            h5path (string):
                Internal path for data you want from a particular hdf5 file.
                Use '{detector}', '{event}' as wildcards to 
                be filled in with the appropriate values for this event. 
                If {detector} is used as a wildcard, results for each appropriate
                detector will be returned.

            datatype (string):
                Can be one of ('array', 'value', 'attribute').
                If it is an attribute, you must also provide attr_name='name-of-attribute'

        Returns:
            A dictionary containing the asked for objects for
            each detector's PSD.

            If detector is not specified:
                A dictionary labelled by detector name, containing objects 
                for each detector

            If detector is specified as a string (e.g. detector='H1'):
                returns the object for that one detector value
        """
        # Run over the list of detectors if not provided
        if ((detectors is None) and ("{detector}" in h5path)):
            detectors = self.SD_ref.available_detectors(self.name)
        
        # Method to interpolate the path
        if replacement_dict is None:
            replacement_dict = {'detector': detectors,  'event': self.name}

        # if the detectors are a list, return a dictionary of results labelled by detector
        if isinstance(detectors, list):
            result = {}
            for ifo in detectors:
                scheme = {'type': datatype, 'name': attr_name, 'path': h5path}
                file = filepath = f"{self.folder}/{event}.hdf5"
                new_replacement_dict = replacement_dict.copy()
                new_replacement_dict['detector'] = ifo
                with h5py.File(file, 'r') as f:
                    result[ifo] = self.SD_ref.read_data_from_file(f, scheme, new_replacement_dict)
        else:
            scheme = {'type': datatype, 'name': attr_name, 'path': h5path}
            file = filepath = f"{self.folder}/{event}.hdf5"
            with h5py.File(file, 'r') as f:
                result = self.SD_ref.read_data_from_file(f, scheme, replacement_dict)
        return result

    def read_posterior_file_from_schema(self, data_name, detectors=None, approximant=None):
        """
        If you've already updated a schema for the posterior database, you can simply just call
        it by name here.

        Example:
        >> df.update_posterior_schema({'calibrations': {'path': '/{approximant}/priors/calibration/{detector}', 
                                                        'type':  'array'}
                                        })
        >> df.read_posterior_file_from_schema('calibrations')
        >> ## Outputs calibration arrays for each detector
        """
        if ((detectors is None) and ("{detector}" in self.PD_ref.schema[data_name]['path'])):
            detectors = self.PD_ref.available_detectors(self.name)

        if approximant is None:
            approximant = self.PD_ref.choose_approximant(self.name)

        if isinstance(detectors, list):
            result = {}
            for ifo in detectors:
                result[ifo] = self.PD_ref.read_data(self.name, data_name, approximant=approximant, detector=ifo)
        else:
            result = self.PD_ref.read_data(self.name, data_name, approximant=approximant, detector=detectors)
        return result


    def read_strain_file_from_schema(self, data_name, detectors=None, approximant=None):
        """
        If you've already updated a schema for the posterior database, you can simply just call
        it by name here. 

        As an example we want to extract the number of time samples in the strain array. This is 
        stored as an attribute names "Npoints" at the path the strain is stored
        Example:
        >> df.update_strain_schema({'Number of points': {'path': '/{detector}/strain/Strain', 
                                                        'type':  'attribute',
                                                        'name': 'Npoints'}
                                        })
        >> df.read_strain_file_from_schema('Number of points')
        >> ## Outputs number of strain time samples for each detector
        """
        if ((detectors is None) and ("{detector}" in self.SD_ref.schema[data_name]['path'])):
            detectors = self.SD_ref.available_detectors(self.name)

        if approximant is None:
            approximant = self.PD_ref.choose_approximant(self.name)

        if isinstance(detectors, list):
            result = {}
            for ifo in detectors:
                result[ifo] = self.SD_ref.read_data(event=self.name, data_name=data_name, detector=ifo)
        else:
            result = self.SD_ref.read_data(event=self.name, data_name=data_name, detector=detectors)
        return result

        


        