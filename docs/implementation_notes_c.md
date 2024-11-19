# 3D SIM algorithm and simrecon.c implementation notes

## definitions

**directions**: Number of illumination pattern directions.  Usually 3.

**phase**: Number of phases for each direction.  Usually 5.

**order**: The separated bands.  Usually 3.

**angle**: The angle of the pattern vector in frequency space. 

**magnitude**:  The length of the pattern vector in frequency space, inverse of the spacing of the pattern

**amplitude**:  (Todo: clarify wording) The modulation depth of the pattern. 

## Summary of code

All processing code is in a single file simrecon.c.

Main simrecon routine starts around line 69.

### File IO

Files are opened using IVE Library starting around line 150 to line 350

Optionally can use .tif reader. To use tiff reader need to set following preprocessor directive 

```
__SIRECON_USE_TIFF__
```

### Threads

Number of threads to use is a parmeter and set at line 370

### Memory allocation and FFT setup

This is done from line 380 to about line 500.

### Main processing loop

begins at about line 522

### drift detection

done at line 550 to about 600

### makematrix

Generates the matrix that is used to separate the raw data into the different bands of the sample information.

### Parameter estimation

Parameter estimation is done as follows

->Loop through each direction
->Loop through each z
loadandflatfield
apodize
rescale
->Loop through each phase
  -> 2D FFT
separate banads

for each direction call

findk0(...) (find initial estimate of modulation wave vector k0 by cross correlation)
fitk0andmodulationamplitudes

### filterbands

Prepares the bands to be assembled:  Applies appropriate OTF-based weighting in the overlaps, uses a Wiener like filter.

- called for each direction

### assemblereadspacebands

Assembles bands in real space.  Calls move.  ```assemblerealspacebands``` is called once for each direction. 

#### move(...)

Moves the contents of the complex fourier space array of size ```(nx/2+1)*ny``` to a bigger array of size ```(zoomfact*nx)*(zoomfact*ny)```.  In most cases ```zoomfact``` is 2. 
- called for each direction

#### After ```move``` is called

- for 3D, each direction has 3 orders.  Order 0 remains at the origin (or frequency space).  Order 1 is shifted by ```mag```.  Order 2 is shifted by 2*```mag```.  