import os
import subprocess
from . import DataFrameClasses as ringdown
from ringdown import PowerSpectrum
from tqdm import tqdm
from .peak import complex_strain_peak_time_td
from . import File
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

approximant_order = ["IMRPhenomPv2",
"IMRPhenomPv3",
"IMRPhenomXPHM",
"IMRPhenomHM",
"IMRPhenom",
"AlignedSpin",
"PrecessingSpin",
"TaylorF2",
"SEOBNRv",
"NRSur7dq4"]

default_schema = {'samples' : {'type': 'array', 'path': '{approximant}/posterior_samples'},
		  'psd': {'type': 'array', 'path': '{approximant}/psds/{detector}'}}


import lalsimulation as ls
import h5py

class PosteriorDatabase:
    def __init__(self, folder, url_df, psd_url_df, strain_url_df, schema=default_schema, approximant_order=approximant_order, cosmo=True):
        self.url_df = url_df
        self._folder = folder
        self.schema = schema
        self.approximant_order = approximant_order
        self.cosmo = cosmo
        self.psd_url_df = psd_url_df
        self.strain_url_df = strain_url_df

    @property
    def folder(self):
        return self._folder

    @folder.setter
    def folder(self, folder):
        if folder[-1] == '/':
            folder = folder[:-1]
        self._folder = folder
        
    def available_detectors(self, event):
        return list(self.strain_url_df[self.strain_url_df.Event == event].Detector.unique())
    
    def event_path(self, event):
        file_type = self.get_file_extension(event)
        post_filename = f"{self.folder}/{event}.{file_type}"
        return post_filename

    def in_catalog(self, event, catalog):
        is_in_the_catalog = self.url_df[self.url_df.event == event].catalog.values[0] == catalog
        return is_in_the_catalog
    
    def in_GWTC1(self, event):
        url = self.url_df[self.url_df.event == event].url.values[0]
        events = self.url_df[self.url_df.url == url].event.unique()
        return ('GW150914' in events)
    
    def in_GWTC3(self, event):
        is_in_the_catalog = self.url_df[self.url_df.event == event].catalog.values[0] == 'GWTC-3'
        return is_in_the_catalog
    
    def event_exists(self, event):
        return os.path.exists(self.event_path(event))
    
    def available_approximants(self, event):
        if self.event_exists(event):
            with h5py.File(self.event_path(event),'r') as f:
                approximants = list(f.keys())
        else:
            print("Hasn't been downloaded")
        approximants2 = [a for a in approximants if a not in ['combined','history', 'version']]
        return approximants2
    
    def choose_approximant(self, event):
        # Chooses the right available waveform approximant based on the priority list provided
        # in self.approximant_order
        list_of_approximants = self.available_approximants(event)
        for test_approx in self.approximant_order:
            options = [approx for approx in list_of_approximants if (test_approx in approx)]
            if len(options) != 0:
                return min(options)
        return None
        
    @property
    def events_present(self):
        files = os.listdir(self.folder)
        files = [f.split('.')[0] for f in files if f[0] != "."] # Keep only those that aren't hidden
        return files
    
    def get_url(self, event):
        mask = (self.url_df.event == event) & (self.url_df.cosmo.isnull() | (self.url_df.cosmo==self.cosmo))
        url = self.url_df.loc[mask,'url'].values[0]
        return url

    def get_file_extension(self, event):
        # This function returns the expected file extension of the finally saved file
        # Determine the file type we need that's inside
        if self.in_catalog(event, 'GWTC-2') or self.in_catalog(event, 'GWTC-2.1'):
            file_type = 'h5'
        elif self.in_catalog(event, 'GWTC-1'):
            file_type = 'h5'
        else:
            file_type = 'h5'
        return file_type

    @staticmethod
    def hdf5_preprocess(df):
        """
        Taken from 
        https://stackoverflow.com/questions/30773073/save-pandas-dataframe-using-h5py-for-interoperabilty-with-other-hdf5-readers

        Convert a pandas DataFrame object to a numpy structured array.

        :param df: the data frame to convert
        :return: a numpy structured array representation of df
        """

        def type_rules(the_type):
            if 'numpy.object_' in str(the_type):
                return np.dtype('S64')
            else:
                return the_type

        v = df.values
        cols = df.columns 
        types = [(k, type_rules(df[k].dtype.type)) for (i, k) in enumerate(cols)]
        dtype = np.dtype(types)
        z = np.zeros(v.shape[0], dtype)
        for (i, k) in enumerate(z.dtype.names):
            z[k] = v[:, i]
        return z, dtype

    
    def download_file(self, event):
        # Get the url to download for this event
        url = self.get_url(event)
        mask = (self.url_df.event == event) & (self.url_df.cosmo.isnull() | (self.url_df.cosmo==self.cosmo))
        file_type = self.url_df.loc[mask,'filename'].values[0].split('.')[-1]
        
        # Check to see if this file is a single event's file or does this compressed file hold multiple events?
        events = self.url_df[self.url_df.url == url].event.unique()
        
        if 'GW150914' in events:
            print("Samples from GWTC-1 come in one file, and there isn't a straightforward")
            print("a way to download only one particular event from GWTC-1.") 
            print("So we have to download them together")
        
        # Download the file
        thefile = File.from_url(url, self.folder, new_filename=f"{event}.{file_type}")
        
        multiple_events = (len(events) > 1)
        
        # Extract and place in the correct location if needed
        if file_type == 'tar':
            # Extract and delete the compressed file
            thefile.extract_here()
            thefile.delete()
            
            # Get the eventual file extension
            file_type = self.get_file_extension(event)

            # Empty out the hdf5 from the extracted folder then delete the folder
            filepath = f"{self.folder}/{event}/{event}.{file_type}"
            new_path = f"{self.folder}/{event}.{file_type}"
            subprocess.run(['mv',filepath, new_path])
            subprocess.run(['rm','-rf',f"{self.folder}/{event}"])
            
        elif file_type == 'zip':
            # Extract and delete the compressed file
            thefile.extract_here()
            thefile.delete()
            
            # GWTC-1 comes in a zip file which always opens up as a folder called pesummary_samples
            # This line is to rename it to a folder with the name of the current event
            if "pesummary_samples" in os.listdir(self.folder):
                subprocess.run(['mv',f"{self.folder}/pesummary_samples",f"{self.folder}/{event}"])
                
            # Empty out the hdf5 from the extracted folder then delete the folder
            downloaded_folder = f"{self.folder}/{event}"
            for file in os.listdir(downloaded_folder):
                file_type = file.split(".")[-1]
                file_eventname = [e for e in events if e in file][0]
                subprocess.run(['mv',f"{downloaded_folder}/{file}", f"{self.folder}/{file_eventname}.{file_type}"])
            subprocess.run(['rm','-rf',f"{self.folder}/{event}"])
            
        # Find the corresponding file in the folder as a sanity check:
        filename = [file for file in os.listdir(self.folder) if event in file][0]
        
        # Return a file with the same 
        return File(f"{self.folder}/{filename}")
    
    def posteriors(self,eventname, peaks=False, f_ref=20.0, f_low=20.0):
        # Download the file if it doesn't exist
        if not self.event_exists(eventname):
            self.download_file(eventname)
            
        replace_names = lambda x: x.replace('C01:','').replace(':HighSpin','').replace('-HS','')
            
        # Create a dataframe of posteriors from the hdf5 or h5 file, or 
        # do the same from the .dat files
        # This adds in the waveform used to generate the posteriors and it's
        # corresponding lalsimulation waveform code
        post_filename = self.event_path(eventname)
        file_type = post_filename.split('.')[-1]
        if (file_type == 'h5') or (file_type == 'hdf5'):
            with h5py.File(post_filename,'r') as f:
                approx = self.choose_approximant(eventname)
                posterior_path = f"/{approx}/posterior_samples"
                all_posteriors = f[posterior_path][:]
                df_posteriors_all = pd.DataFrame(all_posteriors)
                waveform_name = replace_names(approx)
                waveform_code = getattr(ls,waveform_name)
                df_posteriors_all['waveform_name'] = waveform_name
                df_posteriors_all['waveform_code'] = int(waveform_code)
                if f"/{approx}/meta_data/meta_data" in f:
                    f_ref = f[f"/{approx}/meta_data/meta_data/f_ref"][()]
                    f_low = f[f"/{approx}/meta_data/meta_data/f_low"][()]
        elif (file_type == 'dat'):
            df_posteriors_all = pd.read_csv(post_filename,delimiter='\t')
            df_posteriors_all['waveform_name'] = 'IMRPhenomPv2'
            df_posteriors_all['waveform_code'] = int(ls.IMRPhenomPv2)

        # Edit the posteriors so all of them have atleast two columns called
        # 'final_spin' and 'final_mass'. The following just says that if
        # there are non_evolved quantities with no standard counterparts, just
        # create a new column with 'non_evolved' removed
        for x in df_posteriors_all.columns:
            if '_non_evolved' in x:
                new_col = x.replace('_non_evolved','')
                if new_col not in df_posteriors_all.columns:
                    df_posteriors_all[new_col] = df_posteriors_all[x]


        if peaks:
            peak_df = self.calculate_t_peaks(eventname, f_ref=f_ref, f_low=f_low)
            new_df = pd.concat([peak_df, df_posteriors_all],axis=1)
            masses_match = np.all(np.isclose(new_df['final_mass'], new_df['final_mass_check']))
            spins_match = np.all(np.isclose(new_df['final_spin'], new_df['final_spin_check']))
            if masses_match and spins_match:
                return new_df
            else:
                print(masses_match, spins_match)
                print(new_df[['final_mass','final_mass_check']])
                print(new_df[['final_spin','final_spin_check']])
                raise ValueError("The order of the peak times have gone out of sync with the posteriors")

        return df_posteriors_all

    def calculate_t_peaks(self, event, f_low=20.0, f_ref=20.0, recalculate=False):

        # Create the t_peak directory if not already there
        if not os.path.exists(f"{self.folder}/PeakTimes"):
            subprocess.run(["mkdir", f"{self.folder}/PeakTimes"])

        # Check if the file is there, if it is, just read and return
        # it
        if (not recalculate) and os.path.exists(f"{self.folder}/PeakTimes/{event}.csv"):
            return pd.read_csv(f"{self.folder}/PeakTimes/{event}.csv")


        # Get the posterior samples:
        df_posteriors_all = self.posteriors(event)

        print("Calculating the peak times for each sample. May take several minutes.")
        print("Is only calculated once the first time you ask for it for a particular event")
        print(f"Calculating for the event {event}")

        # Calculate all t_peaks
        calculated_times = []

        for i,x in tqdm(df_posteriors_all.iterrows()):
            my_dict = {}
            try:
                t_peak, t_dict, _, _ = complex_strain_peak_time_td(x.to_dict(), wf=int(x['waveform_code']),
                                                                 dt=(1/4096), f_low=f_low, f_ref=f_ref)
                t_dict = {(ifo+"_peak"):v for ifo,v in t_dict.items()}
                my_dict.update(t_dict)
                my_dict.update({"sample_index": i, 'final_mass_check': x['final_mass'], 'final_spin_check': x['final_spin']})
            except:
                t_dict = {(ifo+"_peak"):0.0 for ifo in ['geocent','H1','L1','V1']}
                my_dict.update(t_dict)
                my_dict.update({"sample_index": i, 'final_mass_check': x['final_mass'], 'final_spin_check': x['final_spin']})
                print(f"Predictions failed for event {event}")
            calculated_times.append(my_dict)

        df_times = pd.DataFrame(calculated_times).rename({'geocent_peak': 't_peak'},axis=1)
        df_times.set_index('sample_index',inplace=True)
        df_times.to_csv(f"{self.folder}/PeakTimes/{event}.csv")
        return df_times

    def psd(self, event, detector=None):
        # Return all detectors if none available
        if detector is None:
            detector = self.available_detectors(event)

        # If the event is GWTC-1, there is a seperate PSD file that needs to be downloaded
        if self.in_GWTC1(event):
            url = self.psd_url_df.loc[self.psd_url_df.event == event, 'url'].values[0]
            file_type = url.split('.')[-1]
            filename = f"{event}_psd.{file_type}"
            filepath = f"{self.folder}/{filename}"
            if not os.path.exists(filepath):
                # Download the psd file
                thefile = File.from_url(url, self.folder, new_filename=filename)

            psd_samples = pd.read_csv(filepath, delimiter='\t')
            renaming = {'# Freq (Hz)': 'freq', 'LIGO_Hanford_PSD (1/Hz)': 'H1', 'LIGO_Livingston_PSD (1/Hz)': 'L1', 'Virgo_PSD (1/Hz)': 'V1'}
            if isinstance(detector, list):
                cols = ['freq'] + detector
                psd_samples = psd_samples.rename(renaming, axis=1)
                psd_samples = psd_samples[cols]
                psd_dict = {ifo: PSD(psd_samples[ifo].values, index=psd_samples['freq'].values) for ifo in detector}
                return psd_dict
            else:
                cols = ['freq'] + [detector]
                psd_samples = psd_samples.rename(renaming, axis=1)[cols]
                psd = PSD(psd_samples[detector].values, index=psd_samples['freq'].values)
                return psd
        
        # If the event is not in GWTC-1 then the samples are available in the posterior files.
        # Download the file if it doesn't exist
        if not self.event_exists(event):
            self.download_file(event)
        
        if isinstance(detector, list):
            psd_dict = {}
            file_type = self.get_file_extension(event)
            filepath = f"{self.folder}/{event}.{file_type}"
            for ifo in detector:
                if self.check_data_exists(event, 'psd', detector=ifo):
                    psd_vals = self.read_data(event, 'psd', detector=ifo)
                    psd_dict.update({ifo: PSD(psd_vals[:,1], index=psd_vals[:,0])})
                else:
                    print(f"The PSD for event {event} and detector {ifo} doesn't exist")
            return psd_dict
        else:
            psd_vals = self.read_data(event, 'psd', detector=detector)
            return PSD(psd_vals[:,1], index=psd_vals[:,0])

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

    def check_data_from_file(self, file, scheme, replacement_dict):
        path = self.preprocess_path(scheme['path'], replacement_dict)
        return (path in file)
    
    def read_data(self, event, data_name, approximant=None, detector=None):
        file_type = self.get_file_extension(event)
        file = f"{self.folder}/{event}.{file_type}"
        if approximant is None:
            approximant = self.choose_approximant(event)
        replacement_dict = {'event': event, 'approximant': approximant, 'detector':detector}
        with h5py.File(file, 'r') as f:
            scheme = self.schema[data_name]
            result = self.read_data_from_file(f, scheme, replacement_dict)
        return result

    def check_data_exists(self, event, data_name, approximant=None, detector=None):
        file_type = self.get_file_extension(event)
        file = f"{self.folder}/{event}.{file_type}"
        if approximant is None:
            approximant = self.choose_approximant(event)
        replacement_dict = {'event': event, 'approximant': approximant, 'detector':detector}
        with h5py.File(file, 'r') as f:
            scheme = self.schema[data_name]
            result = self.check_data_from_file(f, scheme, replacement_dict)
        return result




class PSD(PowerSpectrum):
    def __init__(self, *args, **kwargs):
        super(PSD, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return PSD

    def plot(self):
        """
        Do a loglog plot of the powerspectrum
        """
        fig, ax = plt.subplots()
        ax.loglog(self)
        plt.show()

    @property
    def bin_sizes(self):
        """
        Plot the bin sizes of each frequency bin, to detect if there
        are any inconsistencies in the bin sizes
        """
        return pd.Series(self.freq[1:-1] - self.freq[0:-2], index=self.freq[0:-2])

    def low_pad(self, from_freq=None, to_freq=0.0, val=None):
        """
        Pads the power spectrum from from_freq to to_freq
        with the value val or
        """
        # Set parameter values
        df = self.delta_f
        freq = self.freq
        from_freq = from_freq or freq.min()

        # Set fill value and lowest new freq bin
        fill_value = val or self[from_freq]
        new_lowest_freq_bin = to_freq if to_freq < freq.min() else freq.min()

        new_index = np.append(np.arange(new_lowest_freq_bin, from_freq, df), freq)
        new_index = np.array(list(set(new_index)))
        new_index.sort()

        # Reindex the series and replace the value with the set value
        # or the default value, which is the one at freq_from
        a = self.copy()
        a[a.freq < from_freq] = fill_value
        return a.reindex(new_index, fill_value=fill_value)

    def high_pad(self, from_freq=None, to_freq=None, val=None, inclusive=True):
        """
        Pads the power spectrum in such way that all frequency bins 
        above from_freq get assigned the power at from_freq.
        
        The default behaviour is to extend the frequency indices all the 
        way up to the next power of two
        """
        # Set parameter values
        df = self.delta_f
        freq = self.freq
        from_freq = from_freq or freq.max()
        df_end = df if inclusive else 0.0

        # Set fill value and lowest new freq bin
        fill_value = val or self[from_freq]
        next_pow_of_2 = int(2**np.ceil(np.log(freq.max())/np.log(2)))
        new_highest_freq_bin = next_pow_of_2 if to_freq is None else to_freq

        # Create new index
        new_index = np.append(freq, np.arange(from_freq, new_highest_freq_bin + df_end, df))
        new_index = np.array(list(set(new_index).union(set([new_highest_freq_bin]))))
        new_index.sort()

        # Reindex the series and replace the value with the set value
        # or the default value, which is the one at freq_from
        a = self.copy()
        a[from_freq::] = fill_value 
        return a.reindex(new_index, fill_value=fill_value)



