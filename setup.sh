#!/bin/bash

# setup virtual env
virtualenv .env
source .env/bin/activate

# install dependencies
pip install -r requirements.txt
cd src/pydstool
python setup.py install