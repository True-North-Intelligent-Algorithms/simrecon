## Implementation notes Python

The Python wrapper and utilities for processing tiled data are found in the following location

```
.\scripts\Scripts\simrecon_utils.py
```

or see [here](https://github.com/True-North-Intelligent-Algorithms/simrecon/blob/main/scripts/Scripts/simrecon_utils.py)

### Key functions

### simrecon

This is the wrapper to the c code (```sirecon.exe```).  See [here](https://github.com/True-North-Intelligent-Algorithms/simrecon/blob/main/scripts/Scripts/simrecon_utils.py#L730)

Note:  Currently the location of the executable needs is manually set in the file around line 876.  See [here](https://github.com/True-North-Intelligent-Algorithms/simrecon/blob/main/scripts/Scripts/simrecon_utils.py#L876).

### split_process_recombine

This is the function that splits, then processes and then recombines the images.   See [here](https://github.com/True-North-Intelligent-Algorithms/simrecon/blob/main/scripts/Scripts/simrecon_utils.py#L1521)
