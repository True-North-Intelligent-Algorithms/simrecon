## Input parameters  

- **`-ndirs N`**  
  Number of directions in data (default is 3)

- **`-nphases N`**  
  Number of phases in data (default is 5)

- **`-2lenses`**  
  Data acquired with 2 opposing objectives

- **`-bessel`**  
  Data acquired with Bessel beam SIM

- **`-fastSIM`**  
  Data acquired with fast live SIM

- **`-noff N`**  
  Number of switch-off images in nonlinear SIM data (default is 0)

- **`-recalcarray`**  
  How many times to recalculate overlapping arrays (default is 1)

- **`-inputapo N`**  
  Number of pixels to apodize the input data (`-1` to cosine apodize)

- **`-forcemodamp f1 f2... f_norders`**  
  Force the modulation amplitude to be `f1`, `f2`, etc.  
  *Note: If using other than 3 phases, the `-nphases` flag must be specified **before** the `-forcemodamp` flag.*

- **`-nok0search`**  
  Do not search for the best `k0`

- **`-nokz0`**  
  Do not use the `kz0` plane in `makeoverlaps()` or assembly (for bad `kz0`)

- **`-k0searchAll`**  
  Search for `k0` for every time point in a time series

- **`-k0angles f0 f1... f_(ndirs-1)`**  
  User-supplied pattern angle list  
  *Note: The `-ndirs` flag must be specified **before** this flag.*

- **`-k0spacings f0 f1... f_(ndirs-1)`**  
  User-supplied spacings list  
  *Note: The `-ndirs` flag must be specified **before** this flag.*

- **`-fitonephase`**  
  In 2D nonlinear SIM for orders > 1, modulation amplitude's phase will be order 1 phase multiplied by order. Default is using fitted phases for all orders.

- **`-noapodizeout`**  
  Do not apodize the output data

- **`-gammaApo`**  
  Apodize the output data with a power function

- **`-zoomfact`**  
  Factor by which to subdivide pixels laterally

- **`-zzoom`**  
  Factor by which to subdivide pixels axially

- **`-zpadto`**  
  Total sections after zero padding

- **`-explodefact`**  
  Factor by which to exaggerate the order shifts for display

- **`-nofilteroverlaps`**  
  (Used with `-explodefact`) Leave orders unfiltered in overlapping regions (default is to filter)

- **`-nosuppress`**  
  Do not suppress singularity at OTF origins

- **`-suppressR`**  
  Radius of suppression range

- **`-dampenOrder0`**  
  Dampen order 0 in filter bands

- **`-noOrder0`**  
  Do not use order 0 in assembly

- **`-noequalize`**  
  No equalization of input data

- **`-equalizez`**  
  Equalize images across all z-sections and directions

- **`-equalizet`**  
  Equalize all time points based on the first time point

- **`-wiener f`**  
  Set Wiener constant to `f` (default is 0.00)

- **`-wienerInr f`**  
  Set the Wiener constant of the final time point to `wiener` + this value (default is 0.00)

- **`-background f`**  
  Set constant background intensity (default is 515)

- **`-bgInExtHdr`**  
  Background of each section is recorded in the extended header's 3rd float (in Lin's scope)

- **`-otfRA`**  
  Use radially averaged OTFs

- **`-driftfix`**  
  Estimate and fix drift in 3D

- **`-driftHPcutoff x`**  
  Cutoff frequency (fraction of lateral resolution limit) used in high-pass Gaussian filter for drift estimation (default is 0.0)

- **`-fix2Ddrift`**  
  Correct 2D drifts for each exposure within a pattern direction

- **`-fixphasestep`**  
  Correct phase used in separation matrix based on within-direction drift correction

- **`-noff`**  
  Number of switch-off images in nonlinear SIM data

- **`-usecorr corr_file`**  
  Use correction file to perform flatfielding

- **`-nordersout x`**  
  Number of orders to use in reconstruction. Use if different from `(n_phases + 1) / 2`

- **`-angle0 x`**  
  Starting angle (in radians) of the patterns

- **`-negDangle`**  
  Use negative angle step

- **`-ls x`**  
  Line spacing of the illumination pattern (in microns)

- **`-na x`**  
  Effective numerical aperture (NA) of the objective

- **`-nimm x`**  
  Refractive index of the immersion liquid

- **`-saveprefiltered filename`**  
  Save separated bands into a file

- **`-savealignedraw filename`**  
  Save drift-corrected raw images into a file

- **`-saveoverlaps filename`**  
  Save overlaps from `makeoverlaps()` into a file

- **`-nThreads`**  
  Number of threads c-code will use (set to 1 if higher level (ie python with tiles) code will be threaded)

- **`-keeporder2`**  
  If set order 2 coefficients will be used in the order 2 passband instead of mixing with other orders

- **`-help` or `-h`**  
  Print this help message
```
