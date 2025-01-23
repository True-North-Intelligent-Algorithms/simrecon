Build instructions:

1. The "Makefile" in this folder is all you need to build executables on any platform (Linux, Mac OS, Windows). Just cd into this folder and type "make" (or "make clean" first and then "make" if you want to re-build). But before doing that, complete the following step first.

2. For both Linux and Mac OS, you should first download the UCSF Priism package from http://msg.ucsf.edu/IVE/Download/index.html and follow the instruction to install it on your workstation. If all went well, you should have the Shell environment variable, IVE_EXE, defined. If not, go into the Makefile and modify the line

    IVE_BASE     = $(IVE_EXE)/..

to:
    IVE_BASE     = base_folder_name_where_Priism_is_installed

3. For both Linux and Mac OS, you need to have FFTW3 binary and development packages installed. Use the standard package maintenance system on Linux. On Mac, a useful source for these kind of libraries is MacPorts. The Makefile assumes that, on a Mac, all MacPorts libraries are installed under /opt/local/lib, which is the default. If not, modify the line that contains "/opt/local/lib".

4. For Windows 64 bit, install MinGW (http://www.mingw.org) enviroment and mingw-w64 (https://mingw-w64.org/doku.php) compiler toolchain. One can possibly also use Visual Studio, but in this case, I found using MinGW is easier. Had MinGW been an acceptable compiler choice for CUDA, I would have chosen it also for the various CUDA projects.

4.1 When installing MingW, make sure to install "msys-make".
4.2 Install mingw-w64 from https://mingw-w64.org/doku.php/download/mingw-builds . Make sure you select Architecture to be x86-64 at the first installation dialog; also I sometimes have bad luck with "seh" Exception option. If selecting that fails the installation, try change to "sjlj".
4.3 Find the folder "mingw64" inside the folder the last step created (default is C:\Program Files\mingw-w64\x86_64-..., depending on the version number etc.) and copy that whole folder to C:\MinGW
4.4 Use your favorite editor to create a file named ".profile" (note the dot) in your MSys home folder (usually in C:\MingW\msys\1.0\home\your_user_id) and add the following line to it:

PATH=/mingw/mingw64/bin:$PATH

Save it and start Msys from the desktop. From MSys terminal, type "which gcc" and see if it returns with "/mingw/mingw64/bin/gcc.exe". If it does, then proceed.

4.5 Copy the folder "IVE" from this repo to your MSys home folder.

4.6 Either download the TIFF, FFTW3, and LAPACK/BLAS source code and build and install them from inside MSys terminal, or create a folder /usr/local from MSys by typing "mkdir /usr/local", and copy the content of "usr_local" folder to that folder by typing "cp -a usr_local/* /usr/local" in MSys

4.7 Then make as in Step 1. The executable sirecon.exe should run in a Windows CMD window given the following DLLs in the same folder as the executable or in the search path:

libgfortran-3.dll
libgomp-1.dll
libfftw3f-3.dll

They should be all under C:\MinGW folder; look for them there.

Use the following command line in MSys to find out DLL dependencies:
objdump -p DLL_file_name | grep "DLL Name"

Usually it's a safe bet and not too excessive to copy all DLLs under C:\MinGW\mingw64\bin to the same folder as sirecon.exe.

5. This repository also contains source code for making rotaionally averaged OTF files based on PSFs. The executables are called otf2d and radialft (for 2D and 3D SIM, respectively)

5.1. This is how to use otf2d typically:

otf2d -L 3 -H 10 psf_file otf_file

5.2 To use radialft:

radialft psf_file otf_file -fixorigin 3 10 -background 100 -angle 3.74 -leavekz 13 15 3  -na 1.35 -ifixkr 128

-fixorigin: trying to fix the singularity of 3D OTF amplitude around origin. usually "3 10" for 256^2, or "3 20" for 512^2 PSF images.
-background: usually don't need to provide this; program estimate it automatically.
-angle: the SIM pattern angle (in radian) with which the PSF is taken
-na: the approximate detection NA of the objective
-leavekz: used for cleaning the OTFs (i.e. masking out region outside theoretical OTF support regions). The first 2 numbers correspond to the pixel indices where order-1 OTF intersects kz axis. Typically one would run "radialft" without the "-leavekz" flag to obtain a "raw" in order to have better assessment of the PSF/OTF's quality. One of the key signs in raw OTFs to look for is whether signal outside the OTF support is relatively much lower than inside, for all 3 orders. If reasonable then proceeed to a 2nd round of "radialft" with the "-leavekz" flag; one might have to go back and forth several times till the -leavekz numbers are correctly set such that the masking does not affect the OTF support while also not letting too much outside of support in.
-ifixkr: trying to use kz=+/-1 pixels to average out kz=0 pixels, which sometimes can be very noisy at high kx or ky indices. usually use a number that's aboiut 1/4 of the x-y image size.

-- Sometimes the order-1 OTF can look very dim, in which case also try add "-PIshift" flag to the command line.
