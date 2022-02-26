# ringdb
A gravitational event database that queries and locally saves event strain samples,  detector PSDs and posterior samples.

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
