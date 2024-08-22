## Structured Illumination Microscopy Tests

Build instructions

## c code

For windows 64 install MinGW and Msys UCRT

The build starts by running build.sh which in turn calls cmake

Install ffftw, blas and lapack and modify build.sh to point to the correct location of those dependencies on your system. 

Todo:  Add more detail regarding c-build

## Python code

Either run the below pip commands, or ```pip install requirements.txt``` the ```requirement.txt``` file is at the top level of the repository. 

```
pip install ipykernel
pip install numpy==1.24
pip install matplotlib
pip install tifffile
pip install scipy
pip install requests
pip install tqdm
pip install dask
pip install mrc
pip install scikit-image
pip install dphtools
pip install seaborn
pip install tnia-python
```

In your python program add the location of dhputils and ```simrecon_utils.py```.  For example on the main Hess lab machine the path are...  

```
sys.path.insert(1, 'Y:\Cryo_data2\Data Processing Notebooks')
sys.path.insert(1, 'Y:\Cryo_data2\Data Processing Notebooks\Scripts')
```