## 3D SIM algorithm and simrecon.c implementation notes

### definitions

directions:

phase:

order:

angle:

magnitude:

amplitude:

### filterbands

Prepares the bands to be assebled:  Applies appropriate OTF-based weighting in the overlaps, uses a Wiener like filter.

- called for each direction

### assemblereadspacebands

Assembles bands in real space.  Calls move.  ```assemblerealspacebands``` is called once for each direction. 

#### move(...)

Moves the contents of the complex fourier space array of size ```(nx/2+1)*ny``` to a bigger array of size ```(zoomfact*nx)*(zoomfact*ny)```.  In most cases ```zoomfact``` is 2. 
- called for each direction

#### After ```move``` is called

- for 3D, each direction has 3 orders.  Order 0 remains at the origin (or frequency space).  Order 1 is shifted by ```mag```.  Order 2 is shifted by 2*```mag```.  

