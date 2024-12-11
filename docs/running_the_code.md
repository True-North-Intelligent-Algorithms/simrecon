## Setting up a new machine to run code

1.  Install a flavor of Conda (be careful of licensing, may want to look into miniconda or miniforge)

2.  Create an environment from the ```requirement.txt```

$ conda create --name simrecon_environment --file requirements.txt

3.  Get the c executable (```simrecon.exe```) and dependencies

These can be found at the top level of the code repository in the ```bin``` folder.

Or see [here](https://github.com/True-North-Intelligent-Algorithms/simrecon/tree/main/bin)

4.  Point ```simrecon_utils``` to the location of the executable.    

Currently the location of the executable needs is manually set in the file around line 876.  See [here](https://github.com/True-North-Intelligent-Algorithms/simrecon/blob/main/scripts/Scripts/simrecon_utils.py#L876).

5.  Run one of the template notebooks

Run one of the template notebooks for example ```3D SIM reconstruction notebook template.ipynb```

For example see [here](https://github.com/True-North-Intelligent-Algorithms/simrecon/blob/main/notebooks/3D%20SIM%20reconstruction%20notebook%20template.ipynb)

In these notebooks you will need to add some locations to the python path.  

1.  The location where ```dhputils``` is located (this is the ```scripts``` dir at the root of the github repo)
2.  The location where ```simrecon_utils.py``` is located (this is ```scripts\Scripts```)

For example the current code is set up for either the default Janelia machine or Brian Northan's machine.  For a new machine you will have to add corresponding locations. 

```
computer = 'bnort'

import sys

if computer == 'default':
    sys.path.insert(1, 'Y:\Cryo_data2\Data Processing Notebooks')
    sys.path.insert(1, 'Y:\Cryo_data2\Data Processing Notebooks\Scripts')
elif computer == 'bnort':
    sys.path.insert(1, r'C:\Users\bnort\work\Janelia\code\\simrecon\scripts\Scripts')
    sys.path.insert(1, r'C:\Users\bnort\work\Janelia\code\\simrecon\scripts')
else:
    pass
```

6.  Change the OTF Path and data path

You should see a line similar to below.  Change to corresponding location on local machine. 

```home``` is the location of the data (mrc files)

```otf_path``` is the location of the OTFs

```
if computer == 'default': 
    home = r'Y:\Cryo_data2\488nm comparison TILED VS Non-Tiled at different base kwarg settings'
    otf_path = r'Y:\Seyforth\Data For Brian\Cryo-SIM Scope #2 Data (James System)\PSFs (best PSFs and examples of bad ones)\BEAD 2 - NON-AR 1.2W 25ms retake_20240503_170242 BEST PSF!!\computed_OTF_folder'
elif computer == 'bnort':
    home = r'D:\Janelia\gold_standard_set\position 3 nice cells 561nm 4.74um stack_20240530_153702'
    otf_path = r'D:\Janelia\Data 2024-06-03\PSF-OTF used (Davids set of 4 wavelengths)\201909_19-20_best'
```

7. Change the data path