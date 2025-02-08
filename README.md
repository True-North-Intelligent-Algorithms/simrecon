## Structured Illumination Microscopy Experiments

This repo is used to keep track of code that used to test tiled SIM reconstruction on data acquired with a low NA cryogenic 3D structured illumination SIM. 

The goal of our project is to facilitate 3D SIM processing with tiling, and implement and test variations in the reconstruction scheme.  If you are looking for an easy to install (in Python) fast version of simrecon for everyday results check out [cudasirecon](https://github.com/scopetools/cudasirecon).   If you are interested in using our code to perform processing feel free to ask questions on the [Image.sc forum](Image.sc). 

Some of our example data from the cryo-sim can be found [here](https://www.dropbox.com/scl/fo/m47nmkke6kzl4c4b3rm6q/AG2RkOaRT9oRsdhhZcd7QUI?rlkey=e8y48fa1dnl94f2d22fqjcnfz&st=kljvrwnt&dl=0)

This code base is a derivative of Mats Gustafsson & Lin Shao's 3-beam SIM reconstruction software. 

The algorithm as described in [Gustafsson et al (2008) Biophys. J. 94(12): 4957–4970. doi: 10.1529/biophysj.107.120345](https://pmc.ncbi.nlm.nih.gov/articles/PMC2397368/)

We also derive from  a Python wrapper, written by [David Hoffman](https://github.com/david-hoffman) which implements a tiled version of simrecon.

James Seyforth and Brian Northan also made the following contributions

1.  Organized code into one project, added CMake support and added to github repo. 
2.  Better documentation of code organization and build process.
3.  More sophisticated threading to avoid thread resource contention. 
4.  Misc. fixes for some crash scenarios (long file names)
5.  Added a number of notebooks to analyze and troubleshoot results and PSF. 
6.  Added a notebook that uses Napari to visualize results and parameters at each tile location. 
7.  Added constraints to tiled reconstruction to handle occasional invalid results at a tile (this allowed us to experiment with smaller tiles).
8.  Added keeporder2 mode (more emphasis on order2 information) and ran several tests with 'nofilteroverlaps' mode. 
9.  Added new notebooks and scripts for processing. 
10. Established Gold standard image sets and showed a comparison to the original sim recon implementation at Janelia. 
11. Documented in a presentation the experiments and parameters that were tested for purposes of sharing with the community.

Build instructions for Windows

## c code

See [here](docs/build_process.md) for instructions as to how to build c code. 

## Python code

See [here](docs/running_the_code.md) for instructions as to how to run the python code. 

## Notes on Open3DSIM

We also processed data with [Open3DSIM](https://github.com/Cao-ruijie/Open3DSIM).  We made a few modifications so we could run ```Open3DSIM``` with our data.  Most importantly in Open3DSIM the final result has an absolute value operator applied.  In our version we output the data before the ```abs``` operator to better comare with ```simrecon```.  You can find our modified version [here](https://github.com/True-North-Intelligent-Algorithms/Open3DSIM)


