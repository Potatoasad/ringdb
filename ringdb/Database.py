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

    def event(self, eventname):
        return Event(eventname, self)


class Event:
    def __init__(self, eventname, DB_reference):
        self.name = eventname
        self.DB_ref = DB_reference
        self.PD_ref = DB_reference.PosteriorDB
        self.SD_ref = DB_reference.StrainDB

    def posteriors(self):
        return self.PD_ref.posteriors(self.name)

    def psd(self, detector=None):
        return self.PD_ref.psd(self.name, detector=detector)

    def strain(self, detectors=None, duration=32.0):
        return self.SD_ref.strain(self.name, detectors=detectors, duration=duration)

    def read_posterior_file(self, h5path, datatype='array', attr_name=None, detectors=None, approximant=None, replacement_dict=None):
        # Run over the list of detectors if not provided
        if detectors is None:
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
                result[ifo] = self.PD_ref.read_data_from_file(file, scheme, new_replacement_dict)
        else:
            scheme = {'type': datatype, 'name': attr_name, 'path': h5path}
            file = self.PD_ref.event_path(self.name)
            result = self.PD_ref.read_data_from_file(file, scheme, replacement_dict)
        return result

    def read_strain_file(self, h5path, datatype='array', attr_name=None, detectors=None, replacement_dict=None):
        # Run over the list of detectors if not provided
        if detectors is None:
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
                result[ifo] = self.SD_ref.read_data_from_file(file, scheme, new_replacement_dict)
        else:
            scheme = {'type': datatype, 'name': attr_name, 'path': h5path}
            file = filepath = f"{self.folder}/{event}.hdf5"
            result = self.SD_ref.read_data_from_file(file, scheme, replacement_dict)
        return result

        