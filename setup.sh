#!/bin/bash

# setup virtual env
if virtualenv .env ; then
    if source .env/bin/activate ; then

        # install dependencies
        pip install -r requirements.txt
        cd src/pydstool
        python setup.py install

    else
        echo 'activation of virtual enviroment failed'
        echo 'check if source command works'
    fi
else 
    echo 'creation of virtual environment failed. please install requirement:'
    echo 'pip install virtualenv'
fi