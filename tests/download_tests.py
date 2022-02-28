import pytest
import pandas as pd

from ringdb import Database

def create_db():
	db = Database("./Data")
	db.initialize()
	return db


class TestRingdown:

	test_event = "GW190521"
	test_events = ["GW150914",
				   "GW190521", 
				   "GW190426_190642",
				    "GW200220_061928"]
	non_GWTC1_tests = ["GW190521", 
				   	   "GW190426_190642",
				       "GW200220_061928"]

	def test_database_creation(self):
		create_db()

	def test_event_psd(self):
		db = create_db()

		for test_event in self.test_events:
			event = db.event(test_event)
			assert isinstance(event.psd(), dict)

	def test_event_posterior(self):
		db = create_db()

		for test_event in self.test_events:
			event = db.event(test_event)
			assert isinstance(event.posteriors(), pd.DataFrame)

	def test_event_strain(self):
		db = create_db()
		for test_event in self.test_events:
			event = db.event(test_event)
			assert isinstance(event.strain(), dict)

	def test_custom_stain_schema_update(self):
		db = create_db()

		## Only supported for non GWTC-1 events
		for test_event in self.non_GWTC1_tests:
			event = db.event(test_event)
			db.update_strain_schema({'thepoints': {'type': 'attribute', 'name': 'Npoints', 'path': '{detector}/strain/Strain'}})
			event.read_strain_file_from_schema('thepoints')

	def test_custom_posterior_schema_update(self):
		db = create_db()

		## Only supported for non GWTC-1 events
		for test_event in self.non_GWTC1_tests:
			event = db.event(test_event)
			new_schema = {'calibrations': {'type': 'array', 'path': '/{approximant}/priors/calibration/{detector}'}}
			db.update_posterior_schema(new_schema)
			event.read_posterior_file_from_schema('calibrations')

