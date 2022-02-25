__all__ = []

import pandas as pd
import numpy as np
import ringdown
import h5py
import subprocess
import os
import lalsimulation as ls


from .File import *
from .StrainDatabase import *
from .PosteriorDatabase import *
from . import File
from . import StrainDatabase
from . import PosteriorDatabase

folder_post = "./TestingNew/Data/PosteriorData"

from . import metadb

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from . import metadb  # relative-import the *package* containing the templates

posterior_url_path = pkg_resources.open_text(metadb, 'posterior_urls.csv')
# or for a file-like stream:
strain_url_path = pkg_resources.open_text(metadb, 'strain_urls.csv')

#posterior_url_path = os.path.join(os.getcwd(), 'ringdb/metadb/posterior_urls.csv')
#strain_url_path = os.path.join(os.getcwd(), 'ringdb/metadb/strain_urls.csv')

posterior_url_df = pd.read_csv(posterior_url_path)
strain_url_df = pd.read_csv(strain_url_path)

PD = PosteriorDatabase(folder_post, posterior_url_df)
SD = StrainDatabase(folder_post, strain_url_df)
