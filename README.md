# ringdb
A gravitational event database that queries and locally saves event strain samples,  detector PSDs and posterior samples. - Not completely stable yet, fixing issues. 

## Usage:
Define two folders, one to keep the Posteriors and PSDs and one to keep the strain files. Then initialise the `PosteriorDatabase` class and the `StrainDatabase` class
```python
from ringdb import PosteriorDatabase, StrainDatabase, posterior_urls, strain_urls

PD = PosteriorDatabase("/path/to/data/PosteriorDataFolder", posterior_urls)
SD = StrainDatabase("/path/to/data/StrainDataFolder", strain_urls)
```

Then you can call them using the two most important functions. The functions here will download all the requisite files if they need to be downloaded or just read them if the file already exists. 

```python
strain_data = SD.get_strain("GW150914") #  Returns a dictionary of ringdown.Data objects indexed by the detector name
posterior_data = PD.posteriors("GW150914") # Returns a dataframe of all the posterior samples
```
