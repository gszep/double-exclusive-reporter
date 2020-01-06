conda create -n exrep --file requirements_conda.txt
conda activate exrep
cd src/pydstool
python setup.py install
cd ../..