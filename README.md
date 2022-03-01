# ringdb
A gravitational event database that queries and locally saves event strain samples,  detector PSDs and posterior samples.

## Quick Start

```python
from ringdb import Database

# Initialise the database in an empty folder you want your data stored
db = Database("./Data") 
db.initialize()

# Pick an event of your choosing 
event = db.event("GW190521")

# Query it for psds, posteriors and strain data for each detector 
psd = event.psd()	
posteriors = event.posteriors() 
strains = event.strain()
```

## Installation

#### From pip

This package can be downloaded directly from pip. 

```bash
pip install ringdb
```

_(Will error out for most non `x86-64` setups, use Rosetta if on M1)_

#### Virtual Environments

For the latest version, you can use  `make` to run and install the required things in a virtual environment. Replace `VENV_PATH` with any path you're comfortable with making a folder in.  First import once the venv is setup may be slow, but it will be fine after that. 

```bash
git clone https://github.com/Potatoasad/ringdb
cd ringdb
make install VENV_PATH=/path/to/venv
. /path/to/venv/bin/activate
```

## Usage:
### Initialize
Set a folder in your computer somewhere where data will be saved (if the folder isn't there it will be created). Initialize a `Database` object. 
```python
from ringdb import Database

# Sets up the database inside a folder called Data
db = Database("./Data")
db.initialize()
```

Get an object corresponding to a particular event, and query it for psds, posteriors or strain
```python
first_event = db.event("GW150914")
```

### Query

Now you can query a bunch of things from the `Event` object. The most important are:

- `event.posteriors()` -> pandas Dataframe with all the posterior samples
- `event.psd()` -> dictionary of `ringdown.PowerSpectrum` objects (from the [ringdown package](https://github.com/maxisi/ringdown)) for each detector
- `event.strain()` -> dictionary of `ringdown.Data` objects (from the [ringdown package](https://github.com/maxisi/ringdown)) for each detector

#### Get a PSD
```python
first_psd = first_event.psd() # Returns a dictionary labelled by detectors
type(first_psd['H1']) # Ringdown.PowerSpectrum object
```

#### Get posteriors
```python
first_posteriors = first_event.posteriors() # Returns a pandas dataframe
type(first_posteriors) # pandas DataFrame with all the posteriors
```

#### Get the strain files
```python
first_strain = first_event.strain() # Returns a dictionary labelled by detectors
type(first_strain['L1']) # Ringdown.Data object
```

### Custom Queries
_This doesn't work for GWTC-1 events, but will work for the rest_

If you have a particular object you would like to pull from an hdf5 file that is not covered by the above, you can query it directly.

As an example, if I wanted to pull __prior samples for the total mass__. I know for example in GW200316_215756 they will be found in the directory:
`/C01:IMRPhenomXPHM/priors/samples/total_mass` and it needs to read as an `array`.

In general the path seems to follow `/{approximant}/priors/samples/total_mass`. I can then call that for this event:

```python
event = db.event("GW200316_215756")
total_mass_priors = event.read_posterior_file("/{approximant}/priors/samples/total_mass", datatype='array')
```

The database knows how to interpolate things like `{approximant}`, `{detector}` and`{event}`, which correspond to the waveform name, the detectors and the event name.

If you happen to use a particular query a lot you can save it in the database schema, making subsequent accesses easier:

```python
db.update_posterior_schema({'prior_total_mass': {'path': "/{approximant}/priors/samples/total_mass", 'type':'array'}})

total_mass_priors = event.read_posterior_file_from_schema('prior_total_mass')
```



