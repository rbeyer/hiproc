.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

define DOWNLOAD_PYSCRIPT
import urllib.request, sys

urllib.request.urlretrieve(sys.argv[1],sys.argv[2])
endef
export DOWNLOAD_PYSCRIPT


BROWSER := python -c "$$BROWSER_PYSCRIPT"
DOWNLOAD := python -c "$$DOWNLOAD_PYSCRIPT"


help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache
	rm -fr test-resources

lint: ## check style with flake8
	flake8 kalasiris tests

test: test-resources ## run tests quickly with the default Python
	python setup.py test

test-resources: ## Download what we need for testing
	mkdir test-resources
	$(DOWNLOAD) https://hirise-pds.lpl.arizona.edu/PDS/EDR/PSP/ORB_010500_010599/PSP_010502_2090/PSP_010502_2090_RED5_0.IMG test-resources/PSP_010502_2090_RED5_0.img

coverage: ## check code coverage quickly with the default Python
	coverage run --source PyRISE setup.py test
	coverage report -m
# coverage run --source PyRISE setup.py test
# coverage html
# $(BROWSER) htmlcov/index.html
