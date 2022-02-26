__all__ = []

import pandas as pd
import numpy as np
import ringdown
import h5py
import subprocess
import os
import lalsimulation as ls

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources


from .File import *
from .StrainDatabase import *
from .PosteriorDatabase import *
from .Database import *
from . import File
from . import StrainDatabase
from . import PosteriorDatabase
from . import Database

folder_post = "./TestingNew/Data/PosteriorData"

from . import metadb

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from . import metadb 

posterior_url_path = pkg_resources.open_text(metadb, 'posterior_urls.csv')
strain_url_path = pkg_resources.open_text(metadb, 'strain_urls.csv')
psd_url_path = pkg_resources.open_text(metadb, 'psd_urls.csv')

posterior_urls = pd.read_csv(posterior_url_path)
strain_urls = pd.read_csv(strain_url_path)
psd_urls = pd.read_csv(psd_url_path)

PD = PosteriorDatabase(folder_post, posterior_urls, psd_urls, strain_urls)
SD = StrainDatabase(folder_post, strain_urls)
