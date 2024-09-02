## 3D SIM algorithm and simrecon.c implementation notes

### assemblereadspacebands

Assembles bands in real space.  Call move

#### move(...)

Moves the contents of the complex fourier space array of size ```(nx/2+1)*ny``` to a bigger array of size ```(zoomfact*nx)*(zoomfact*ny)```.  In most cases ```zoomfact``` is 2. 



