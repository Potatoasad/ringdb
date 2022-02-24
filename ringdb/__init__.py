__all__ = []

from .File import *
from .StrainDatabase import *
from .PosteriorDatabase import *
from . import File
from . import StrainDatabase
from . import PosteriorDatabase

folder_post = "./TestingNew/Data/PosteriorData"

posterior_url_df = pd.read_csv('./metadb/posterior_urls.csv')
event_df = pd.read_csv('./metadb/strain_urls.csv')

PD = PosteriorDatabase(folder_post, posterior_url_df)
SD = StrainDatabase(folder_post, event_df)