import os
import subprocess
import ringdown
import pandas as pd
import numpy as np
import h5py
from . import File

default_schema = {'sample' : {'type': 'array', 'path': '{detector}/strain/Strain'},
                  't0': {'type': 'attribute', 'name': 'Xstart', 'path': '{detector}/strain/Strain'},
                  'dt': {'type': 'attribute', 'name': 'Xspacing', 'path': '{detector}/strain/Strain'},
                  'Npoints': {'type': 'attribute', 'name': 'Npoints', 'path': '{detector}/strain/Strain'}
                 }
        
class StrainDatabase:
    def __init__(self, folder, url_df, schema=default_schema):
        self.url_df = url_df
        if folder[-1] == '/':
            folder = folder[:-1]
        self.folder = folder
        self.schema = schema
        
    def available_detectors(self, event):
        return list(self.url_df[self.url_df.Event == event].Detector.unique())
        
    @property
    def events_present(self):
        files = os.listdir(self.folder)
        files = [f.split('.')[0] for f in files if f[0] != "."] # Keep only those that aren't hidden
        return files
    
    def get_url(self, event, detector, duration=32.0):
        mask = (self.url_df.Event == event) & (self.url_df.Detector == detector) & (self.url_df.Duration == duration)
        url = self.url_df.loc[mask,'Url'].values[0]
        return url
    
    def download_file(self, event, detector, duration=32.0):
        url = self.get_url(event, detector, duration)
        thefile = File.from_url(url, self.folder, new_filename=f"{event}-{detector}.hdf5") 
        return thefile
        
    def combine_detector_files(self, event, detector_files):
        # Create a new hdf5 file called {event}.hdf5
        with h5py.File(f'{self.folder}/{event}.hdf5','w') as f:
            for ifo, detector_file in detector_files.items():
                
                # Open the downloaded detector strains
                file = h5py.File(detector_file.path,'r') 
                
                # Copy each file into the newly created hdf5 under an internal path like /H1 or /L1
                h5py.h5o.copy(file.id, b"/", f.id, f"/{ifo}".encode()) 
        
    def make_event_file(self, event, duration=32.0):
        # Download all the events available
        files = {}
        for ifo in self.available_detectors(event):
            files[ifo] = self.download_file(event, ifo, duration=duration)
            
        # Combine them into one file named {event}.hdf5
        self.combine_detector_files(event, files)
        
        # Delete the downloaded detector files:
        for ifo, file in files.items():
            file.delete()
    
    @staticmethod
    def preprocess_path(path, replacement_dict):
        rep = lambda x: '' if x is None else x
        replacements = { ("{"+key+"}"): rep(value) for key,value in replacement_dict.items() }
        for key, value in replacements.items():
            path = path.replace(key,value)
        return path
            
    def read_data_from_file(self, file, scheme, replacement_dict):
        path = self.preprocess_path(scheme['path'], replacement_dict)
        if scheme['type'] == 'attribute':
            l = file[path].attrs[scheme['name']]
        elif scheme['type'] == 'array':
            l = file[path][:]
        elif scheme['type'] == 'value':
            l = file[path][()]
        else:
            l = None
        return l
    
    def read_data(self, event, detector, data_name):
        file = f"{self.folder}/{event}.hdf5"
        replacement_dict = {'event': event, 'detector': detector}
        with h5py.File(file, 'r') as f:
            scheme = self.schema[data_name]
            result = self.read_data_from_file(f, scheme, replacement_dict)
        return result
            
    def strain(self, event, detectors=None, duration=32.0):
        # Download the file if the file doesn't exist
        if event not in self.events_present:
            self.make_event_file(event, duration=duration)
            
        # Prepare which detectors which need to be returned
        if detectors is None:
            detectors = self.available_detectors(event)
            
        # Grab the data you need from the detectors you need
        filepath = f"{self.folder}/{event}.hdf5"
        
        if isinstance(detectors,list):
            # If you pass a list of detectors, or None, you'll get a dictionary of Data objects
            strain = {}
            for ifo in detectors:
                h  = self.read_data(event, ifo, 'sample')
                t0 = self.read_data(event, ifo, 't0')
                dt = self.read_data(event, ifo, 'dt')
                strain[ifo] = ringdown.Data(h, index=t0 + dt*np.arange(len(h)), ifo=ifo)
        else:
            # If you pass just one detector string you'll get one data object
            h  = self.read_data(event, detectors, 'sample')
            t0 = self.read_data(event, detectors, 't0')
            dt = self.read_data(event, detectors, 'dt')
            strain = ringdown.Data(h, index=t0 + dt*np.arange(len(h)), ifo=detectors)
        return strain

