SHELL = /bin/bash

venv: 
	python3 -m venv $(VENV_PATH)
	cd $(VENV_PATH)/bin/
	. activate
	python -m pip install ringdb

.PHONY: install

install:
	{ \
	python3 -m venv $(VENV_PATH) ;\
	. $(VENV_PATH)/bin/activate ;\
	pip install --upgrade pip ;\
	pip install . ;\
	}


sample_test: ./tests/download_tests.py
	pytest ./tests/download_tests.py

full_test: ./tests/full_test.py
	python3 -i ./tests/full_test.py
	
