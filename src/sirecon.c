/***************************** sirecon.c  *********************************/
/*     Uses three 3D image stacks acquired with stripe illumination, with */
/*     stripe phase of nominally 0, 120 and 240 degrees (relative to      */
/*     some unknown phase offset), to generate an image with enhanced     */
/*     resolution in the stripe direction.  Using three such image        */
/*     triples in different directions, nearly isotropic lateral          */
/*     resolution enhancement can be achieved.                            */
/**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include <math.h>

#include <float.h>
#include <string.h>

#ifdef __SIRECON_USE_TIFF__

#define cimg_use_tiff
#include <CImg.h>
using namespace cimg_library;
#else
#include <IMInclude.h>  // MRC file I/O routines
#endif

#include <omp.h>   // OpenMP library for multithreading via preprocessor

#include "sirecon.h"

#ifdef _WIN32
#include <windows.h> // for GetSystemInfo()
#else
#include <unistd.h> // for sysconf()
#endif


// for IMLIB with vs >2015
// https://stackoverflow.com/questions/30412951/unresolved-external-symbol-imp-fprintf-and-imp-iob-func-sdl2
#include <stdio.h>

static FILE* temp[3];
static pthread_mutex_t temp_mutex = PTHREAD_MUTEX_INITIALIZER;

void init_iob() {
    pthread_mutex_lock(&temp_mutex);
    temp[0] = stdin;
    temp[1] = stdout;
    temp[2] = stderr;
    pthread_mutex_unlock(&temp_mutex);
}

FILE** __imp___iob_func() {
    FILE** temp_copy;
    pthread_mutex_lock(&temp_mutex);
    temp_copy = temp;
    pthread_mutex_unlock(&temp_mutex);
    return temp_copy;
}

int nthreads;  /* for multithreaded FFTW */

float BesselNA = 0.672; //.5166; // 
float BesselLambdaEx = .488;

float * boundaryAngles=NULL;


int main(int argc, char **argv)
{
  char ifiles[200], ofiles[200], otffiles[200];
  int nx, ny, nz, nz0, zoffset, nxy, z, norders, nordersIn, order;
  int nphases, ndirs, sizeOTF;
  float *buffer;
  fftwf_complex *bigbuffer=0;
  int i, j, iw, it;
#ifdef __SIRECON_USE_TIFF__
  CImg<> input_tiff, output_tiff, otf_tiff;
#else
  int istream_no=1, ostream_no=2, otfstream_no=3;
  int aligned_stream_no = 10, separated_stream_no=11, overlaps_stream_no=12;
  short nwaves, wave[5], ntimes;
  IW_MRC_HEADER header;
  IW_MRC_HEADER aligned_header;  //used for saving drift-corrected off+exp raw images
  IW_MRC_HEADER sep_header;  // used for saving separated bands
  IW_MRC_HEADER overlaps_header;  // used for saving overlaps
  char titles[1000];
#endif
  float **rawImages, **offImages=NULL, *oneOffImage=NULL, *refImage=NULL;
  float *outbuffer=0;
  fftwf_complex **otf, *overlap0=NULL, *overlap1=NULL;

  fftwf_complex **bands, ***savedbands;
  /* we want savedbands: ndirs pointers to nphases pointers to complex arrays of size nz*ny*(nx/2+1)   */
  float ***bufOffImages;  // like savedbands for switch-off images in NL SIM

  fftwf_complex  **amp, amp_inv, amp_combo;
  float  corr_coeff;
  float *floatsection;
  float dxy, dz, dkx, dky, inscale, *background, backgroundExtra=0, *slope;
  float dist, delta_angle;
  float k0magguess;   // in 1/micron, not in pixels!!
  vector *k0, *k0guess, *k0_time0, deltak0; // in pixels
  vector3d *driftlist=0 /* to be used for section-by-section drift correction, on top of between-dirs drift fix */;
  float *phaseList = 0;
  vector *driftAcq=0; // To store drifts estimated in the acquisition software;
  vector3d *drifts = 0;
  float *phaseAbs=0; // To store absolute phase applied for each exposure in the acquisition software;
  float *timestamps=0;  // acquisition time stamp of each exposure to be read from the extended header
  float *timestamps_for_fitting = 0;
  float *expDose=0; // exposure dose reported in extended header
  vector3d *drift_bt_dirs=0;
  int direction;
  float ampfact;
  float *sepMatrix, *noiseVarFactors;
  int ixyz[3], mxyz[3], pixeltype;      /* variables for IMRdHdr call */
  float min, max, mean;      /* variables for IMRdHdr call */
  int ifilein, ofilein, otffilein;
  fftwf_plan fftplan=NULL, rfftplan, fftplan1D=0;
  /* int skips, compen_beadsize=0; */
  double * sum_dir0_phase0, *sum_dir0_phase0_off;    /* for storing the image sum of phase 0 of dir 0; used in rescale() */

  ReconParams myParams;
  float maxval = -FLT_MAX, minval = FLT_MAX; /* used in output file header writing */

  float *est_wiener;


  /*  Set default values to various parameters before parsing command line arguments */
  SetDefaultParams(&myParams);

  /* Flags to indicate whether in, out, and OTF files are specified in the command line */
  ifilein = 0;
  ofilein = 0;
  otffilein = 0;

  /* Parse the command line, results fill in the fields of myParams */
  if (!commandline(argc, argv, &myParams, ifiles, ofiles, otffiles, &ifilein, &ofilein, &otffilein))
    return -1;

  /* These two variables are very often used; so make a copy of them from myParams: */
  nphases = myParams.nphases;
  ndirs   = myParams.ndirs;
  printf("nphases=%d, ndirs=%d\n", nphases, ndirs);

#ifndef __SIRECON_USE_TIFF__
  /* Suppress IVE display of file headers */
  IMAlPrt(0);
#endif
  /* Open the input raw data file */
  if (!ifilein) {
    printf(" Data filename=");
    fgets(ifiles, 200, stdin);
    ifiles[strlen(ifiles)-1] = '\0';
  }
#ifdef __SIRECON_USE_TIFF__
  input_tiff = CImg<>(ifiles);
#else
  if (IMOpen(istream_no, ifiles, "old")) {
    fprintf(stderr, "File %s does not exist.", ifiles);
    return -1;
  }
#endif

  /* Create output file */
  if (!ofilein) {
    printf(" Output filename=");
    fgets(ofiles, 200, stdin);
    ofiles[strlen(ofiles)-1] = '\0';
  }
#ifdef __SIRECON_USE_TIFF__
  output_tiff = CImg<>(); /* will assign dimensions and all later */
#else
  if (IMOpen(ostream_no, ofiles, "new")) {
#endif
    fprintf(stderr, "File %s can not be created.", ofiles);
    return -1;
  }

  /* Open OTF file */
  if (!otffilein) {
    printf(" OTF filename=");
    fgets(otffiles, 200, stdin);
    otffiles[strlen(otffiles)-1] = '\0';
  }
#ifdef __SIRECON_USE_TIFF__
  otf_tiff = CImg<>(otffiles);
#else
  if (IMOpen(otfstream_no, otffiles, "ro")) {
#endif
    fprintf(stderr, "File %s does not exist.", otffiles);
    return -1;
  }

  /* Retrieve header information of the raw data file */
#ifdef __SIRECON_USE_TIFF__
  nz = input_tiff.depth();
  nx = input_tiff.width();
  ny = input_tiff.height();
  nwaves = input_tiff.spectrum();
  nz /= nwaves * nphases * ndirs * (1+myParams.nOffImages);
  /* Presumably ntimes, wavelengths, dy, dz are passed in from commandline? */
  ntimes = myParams.ntimes;
  dxy = myParams.dy;
  dz = myParams.dz;
  for (iw=0; iw<nwaves; iw++)
    wave[iw] = myParams.wavelengths[iw];
#else
  IMRdHdr(istream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(istream_no, &header);
  nx=header.nx; ny=header.ny;
  nz=header.nz/(header.num_waves*header.num_times);    /* header.nz is defined as the total number of sections = nz*nwaves*ntimes (* ndirs*nphases in our case) */
  nz /= nphases*ndirs*(1+myParams.nOffImages);
  nwaves=header.num_waves;
  ntimes=header.num_times;
  wave[0]=header.iwav1; wave[1]=header.iwav2; wave[2]=header.iwav3;
  wave[3]=header.iwav4; wave[4]=header.iwav5;
  /* dxy: lateral pixel size; dz: axial pixel size; both in microns */
  dxy = header.ylen;
  dz = header.zlen;

  /* Initialize headers for intermediate output files if requested */
  if (myParams.bSaveAlignedRaw)
    memcpy(&aligned_header, &header, sizeof(header));
  if (myParams.bSaveSeparated)
    memcpy(&sep_header, &header, sizeof(header));
  if (myParams.bSaveOverlaps)
    memcpy(&overlaps_header, &header, sizeof(header));
#endif

  if (nz > 1 && myParams.b2D) {
    /* nz>1 && b2D means dataset's each time point has multiple z slices of 2D SIM images */
    /* To handle this and avoid any downstream modifications, we modify the input file's header
       to convert multiple z slices to extra time points and keep nz at 1 */
    IMAlZWT(istream_no, header.nz/(header.num_times * nz), nwaves, header.num_times * nz, ZTW_SEQUENCE);
    /* IMWrHdr(istream_no, header.label, 1, header.amin, header.amax, header.amean); */
    ntimes *= nz;
    nz=1;
  }
 
  if (myParams.nzPadTo)
    nz0 = myParams.nzPadTo;
  else
    nz0 = nz;
  printf("nx=%d, ny=%d, nz=%d, nwaves=%d, ntimes=%d\n", nx, ny, nz, nwaves, ntimes);


  /* Dealing with drift-caused phase step deviation, mainly for 2D nonlinear SIM */
  /* Extended header contains experiment-time drift correction applied, time stamps, etc., that will be used in fixing the phase step */
#ifndef __SIRECON_USE_TIFF__
  if (myParams.bFixPhaseStep) {
    /*
       Up till now, header.nz still holds the number of total exposures in the raw data.
       So the following is ok.
    */
    if (nz * nphases *ndirs *(1+myParams.nOffImages) * sizeof(myExtHeader) > header.inbsym) {
      fprintf(stderr, "The extended header does not seem to have info needed for option -fixphistep\n");

/*       if (nz==1) */
/*         myParams.bFixPhaseStep = 0; */

      // then we need to fake some timestamps
      timestamps = malloc(header.nz * sizeof(float));
      for (i=0; i<header.nz; i++) {
        timestamps[i] = i;
      }
    }
    else {
      // load from extended header acquistion-time phase correction values, time stamps, exposure dose, etc.
      float *ptr;
      void * extHeader;
      driftAcq   = malloc(header.nz * sizeof(vector));
      phaseAbs   = malloc(header.nz * sizeof(float));
      timestamps = malloc(header.nz * sizeof(float));
      expDose    = malloc(header.nz * sizeof(float));

      extHeader = malloc(header.inbsym);
      IMGetExHdr(istream_no, extHeader);
      ptr = (float*) extHeader;
      for (i=0; i<header.nz; i++) {
        /*
          float fields in Hesper's extended header are currently defined as:
          (timestamps, phase corr, exposure dose, drift_x, drift_y).
        */
        timestamps[i] = *ptr++;
        phaseAbs[i]   = *ptr++;
        phaseAbs[i]  *= M_PI / 180.0;   // unit used in the extended header is degrees
        expDose[i]    = *ptr++;
        driftAcq[i].x = *ptr++;
        driftAcq[i].y = *ptr++;
      }
      free(extHeader);
    }
  }
#endif
  /* allocate drift_bt_dirs */
  drift_bt_dirs = (vector3d *) malloc(ndirs * sizeof(vector3d));
  drift_bt_dirs[0].x = 0;
  drift_bt_dirs[0].y = 0;
  drift_bt_dirs[0].z = 0;

  /* Assign norders here so that otf can be allocated correctly */
  nordersIn = nphases/2+1;  /* number of orders of bands (including the center, zeroth, order) */
  /* Sometimes number of output orders is specified and could be less than nordesIn */
  if (myParams.norders_output != 0)
    norders = myParams.norders_output;
  else
    norders = nordersIn;

#ifdef __SIRECON_USE_TIFF__
  determine_otf_dimensions(otf_tiff, norders, nz, &myParams, &sizeOTF);
#else
  determine_otf_dimensions(otfstream_no, norders, nz, &myParams, &sizeOTF);
#endif

  /* allocate space for OTFs */
  otf = (fftwf_complex **) malloc(norders * sizeof(fftwf_complex *));
  for (i=0; i<norders; i++)
    otf[i] =  (fftwf_complex *) malloc(sizeOTF * sizeof(fftwf_complex));

#ifdef __SIRECON_USE_TIFF__
  loadOTFs(otf_tiff, otf, norders, nz, sizeOTF, &myParams);  //helpers.c
#else
  loadOTFs(otfstream_no, otf, norders, nz, sizeOTF, &myParams);  //helpers.c
#endif


  /* Now update the output file's header info; "header" is reused from the data file's header */
  header.mode = IW_FLOAT;
  header.nz = nz * nwaves * ntimes * myParams.z_zoom;
  header.nx *= myParams.zoomfact;
  header.ny *= myParams.zoomfact;
  header.xlen /= myParams.zoomfact;
  header.ylen /= myParams.zoomfact;
  header.zlen /= myParams.z_zoom;
  header.inbsym = 0;  /* Discard all extended header info. Is this the right way? */
  IMPutHdr(ostream_no, &header);
  IMAlCon(ostream_no, 0);

  /* Read in background and CCD correction file */
  background = (float *) malloc(nx*ny*sizeof(float));
  slope = (float *) malloc(nx*ny*sizeof(float));
  if (myParams.bUsecorr) {  // flatfield correction of measured data using calibration data
    printf("loading CCD calibration file\n");
    getbg_and_slope(myParams.corrfiles, background, slope, nx, ny);
  }
  else
    for (i=0; i<nx*ny; i++) {
      background[i] = myParams.constbkgd;   /* use the constant background value given by user */
      slope[i] = 1.0;
    }
  backgroundExtra = 0;


#ifndef NDEBUG
  nthreads = 1;
#else
#ifdef _WIN32
  SYSTEM_INFO siSysInfo;
  GetSystemInfo(&siSysInfo);
  nthreads = siSysInfo.dwNumberOfProcessors;
#else
  nthreads = sysconf(_SC_NPROCESSORS_CONF); /* Get number of processors */
#endif
#endif
  nthreads = 1;
  omp_set_num_threads(nthreads);
  printf("Number of threads used: %d\n", nthreads);


  /* If we are going to use "fftw_complex **bands", we need to allocate (fftw_complex **) memory to bands or rawImages here too */
  savedbands = (fftwf_complex ***) malloc(ndirs * sizeof(fftwf_complex **));
  bufOffImages = (float ***) malloc(ndirs * sizeof(float **));

  for( direction=0; direction<ndirs; direction++) {
    savedbands[direction] = (fftwf_complex **) malloc(nphases*sizeof(fftwf_complex *));
    bufOffImages[direction] = (float **) malloc(nphases*sizeof(float *));
    for (i=0; i<nphases; i++) {
      savedbands[direction][i] = (fftwf_complex *) fftwf_malloc((nx/2+1)*ny*nz0 * sizeof(fftwf_complex));  //nz0 could be >nz if padding
      if (!savedbands[direction][i]) {
        printf("Not enough memory for bands after direction %d\n", direction);
        return -1;
      }

      if (myParams.nOffImages)
        bufOffImages[direction][i] = (float *) malloc((nx+2)*ny*nz0 * sizeof(float));
      else
        bufOffImages[direction][i] = 0;
    }
  }

  nxy = (nx+2) * ny;

  /* FFT plans:
     rfftplan  -- 2D r2c FFT of each input z section (whether 2D or 3D data);
     fftplan1D -- if 3D input, 1D forward c2c FFT through the 2D FFT sections along kz, for each (kx,ky) pixel;
     fftplan   -- backward c2c FFT to obtain the real-space output bands before they get assembled into final output
  */
  if (!fftwf_init_threads()) { /* one-time initialization required to use threads */
    printf("Error returned by fftwf_init_threads()\n");
  }

  printf("Creating FFT plans...\n");
  fftwf_plan_with_nthreads(nthreads);

  rfftplan = fftwf_plan_dft_r2c_2d(ny, nx, (float*) savedbands[0][0], savedbands[0][0], FFTW_ESTIMATE);
  if (nz0>1) // 1D FFT of a stack of 2D FFTs
    fftplan1D = fftwf_plan_many_dft(1, &nz0, nxy/2, savedbands[0][0], NULL, nxy/2, 1,
                                    savedbands[0][0], NULL, nxy/2, 1, FFTW_FORWARD, FFTW_ESTIMATE);

  dkx= 1 / (nx*dxy);
  dky= 1 / (ny*dxy);
  k0magguess = 1.0 / myParams.linespacing;  // in 1/micron

  /* In 2D and 3D SIM, orders mean different things -- order 1 in 2D is actually order 2 in 3D in terms of
     k vector magnitude, for example. To unify the meaning of order 1 in 2D and 3D cases, k0magguess is
     halved in 3D case.
  */
  if (nz>1)
    k0magguess /= nordersIn -1; // convert k0mag into the k0mag corresponding to the lowest excitation harmonic

  /* Initialization of pattern wave vectors*/
  k0 = (vector *) malloc(ndirs * sizeof(vector));
  k0_time0 = (vector *) malloc(ndirs * sizeof(vector));
  k0guess = (vector *) malloc(ndirs * sizeof(vector));

  delta_angle = M_PI/ndirs;
  if (myParams.bNegAngleStep)
    delta_angle *= -1;

  for (i=0; i<ndirs; i++) {
    float k0angleguess;
    if (myParams.k0angles == 0)
      k0angleguess = myParams.k0startangle + i*delta_angle;
    else
      k0angleguess = myParams.k0angles[i];

    // k0guess is in 1/micron
    k0guess[i].x = k0magguess * cos(k0angleguess);
    k0guess[i].y = k0magguess * sin(k0angleguess);
  }


  sepMatrix = (float *) malloc(nphases*(norders*2-1)*sizeof(float));
  noiseVarFactors  = (float *) malloc(ndirs*norders*sizeof(float));

  /* Allocate buffer to load single section of raw indata;
     IVE library takes care of converting uint16 to float. */
  buffer =  (float *) malloc((nx*ny) * sizeof(float));

  for (i=0; i<ndirs*norders; i++)
    noiseVarFactors[i] = 1.0;

  /* The standard evenly-spaced-phase-step separation matrix; can be overwritten later */
  makematrix(nphases, norders, 0, 0, sepMatrix, noiseVarFactors);

  inscale = 1.0/(nx*ny*nz0*myParams.zoomfact*myParams.zoomfact*myParams.z_zoom*ndirs);  /* not sure if 'ndirs' is needed here */

  amp = (fftwf_complex **) malloc(ndirs*sizeof(fftwf_complex*));
  for (direction=0; direction <ndirs; direction ++) {   /* initialize amp[dir] for order 0. */
    amp[direction] = (fftwf_complex *) malloc(norders*sizeof(fftwf_complex));
    amp[direction][0] = 1;
  }

  sum_dir0_phase0 = (double *) malloc(nz*nwaves*sizeof(double));
  sum_dir0_phase0_off = (double *) malloc(nz*nwaves*sizeof(double));

  if (myParams.nOffImages)
    oneOffImage = malloc(nxy*sizeof(float));

  if (myParams.bSaveAlignedRaw) {
    IMOpen(aligned_stream_no, myParams.fileRawAligned, "new");
    if (myParams.nOffImages) {
      aligned_header.num_waves = 2;
      aligned_header.interleaved = WZT_SEQUENCE;
    }
    aligned_header.mode = IW_FLOAT;
    aligned_header.inbsym = 0;
    IMPutHdr(aligned_stream_no, &aligned_header);
  }

  if (myParams.bSaveSeparated) {
    IMOpen(separated_stream_no, myParams.fileSeparated, "new");
    sep_header.nx = (nx+2)/2;    // saved will be separated FFTs
    if (myParams.nOffImages) {
      sep_header.num_waves = 2;
      sep_header.interleaved = WZT_SEQUENCE;
    }
    sep_header.mode = IW_COMPLEX;
    sep_header.inbsym = 0;
    IMPutHdr(separated_stream_no, &sep_header);
  }

  if (myParams.bSaveOverlaps) {
    IMOpen(overlaps_stream_no, myParams.fileOverlaps, "new");
    overlaps_header.nz = nz*2*ndirs*ntimes*nwaves;
    overlaps_header.num_waves = 2;  // save overlap 0 and 1 as wave 0 and 1 respectively
    overlaps_header.interleaved = WZT_SEQUENCE;
    overlaps_header.mode = IW_COMPLEX;   // saved will be full-complex overlaps in real space
    overlaps_header.inbsym = 0;
    IMPutHdr(overlaps_stream_no, &overlaps_header);
  }


  if (myParams.nzPadTo)
    zoffset = (nz0-nz)/2;
  else
    zoffset = 0;


  est_wiener = malloc(ndirs*sizeof(float));

  for (it=0; it<ntimes; it++) {
    int phase;
    int nDirPhases = ndirs * nphases;

    for (iw=0; iw<nwaves; iw++) {

      float rdistcutoff;

      if (it==0 && iw==0) {
        overlap0 = (fftwf_complex *) fftwf_malloc( nx*ny*nz0* sizeof(fftwf_complex) );
        overlap1 = (fftwf_complex *) fftwf_malloc( nx*ny*nz0* sizeof(fftwf_complex) );
        if(!overlap0 || !overlap1) {
          printf("Not enough memory for overlap arrays\n");
          return -1;
        }
      }

      rdistcutoff = myParams.na*2 / (wave[iw]/1000.0);    /* OTF support radial limit in 1/micron */
      printf("rdistcutoff = %f 1/micron\n", rdistcutoff);

      myParams.wiener += myParams.wienerInr/ntimes;
      /* printf("wiener = %f\n", myParams.wiener); */
      
      /* Determine amount of lateral and axial drift */
      if (myParams.bFixdrift || myParams.bFix2Ddrift) {

        // determines the drift between dirs, without using it yet
        if (!myParams.bFastSIM)
          determinedrift_bt_dirs(istream_no, iw, it, drift_bt_dirs, nx, ny, nz, ndirs, nphases,
                                 background, dxy, dz, rdistcutoff, wave[iw], &myParams);
          
        /* The reason to estimate drift here instead of later, where fixdrift_bt_dirs() is called, is that
           only the middle sections of the 3D stack are used in estimating, and the bands will have been 
           already FFT'ed when calling fixdrift_bt_dirs(). */

        for (direction=1; direction<ndirs; direction++)
          printf("Drift for direction %d is (%f, %f, %f) pixels\n", direction,
                 drift_bt_dirs[direction].x, drift_bt_dirs[direction].y, drift_bt_dirs[direction].z);


        if (it==0 && iw==0) {
          driftlist = malloc(nz*nphases*ndirs*sizeof(vector3d));
          drifts = malloc(nphases*sizeof(vector3d));
          phaseList = malloc(nphases*sizeof(float));
        }

        if (nz>1 && myParams.bFix2Ddrift) {
          if (myParams.bFastSIM) {
            /* determinedrift_dirs_interleaved(istream_no, iw, it, driftlist, nx, ny, nz, ndirs, nphases, dy, dz, &myParams); */
          }
          else {
            if (it==0 && iw==0)
              timestamps_for_fitting = malloc(ndirs * sizeof(float));

            /* for 3D, use the middle section of each direction as the time points for fitting */
            for (direction=0; direction<ndirs; direction++)
              timestamps_for_fitting[direction] = timestamps[direction*nz*nphases + nz*nphases/2];

            fitXYdrift(drift_bt_dirs, timestamps_for_fitting, ndirs, driftlist, timestamps, ndirs*nz*nphases);
            /*
              output x-y drift estimates, driftlist, for every section, either by fitting between-dirs drift estimates
              (as is done now) or section-to-section drift estimates. x-y drift estimates will be converted to phase
              later in makematrix(), which also needs acquisition-time phase corrections, which are stored in extended
              header
            */

            /* Subtract the fitted drift by the mean drift of each direction,
               because later fixdrift_bt_dir() will be called whether 2D or 3D */
            for (direction = 1; direction < ndirs; direction ++) {
              int offset = direction*nz*nphases;
              for (z=0; z<nz*nphases; z++) {
                driftlist[offset+z].x -= drift_bt_dirs[direction].x;
                driftlist[offset+z].y -= drift_bt_dirs[direction].y;
              }
            }
          }
        }


        if (it==0 && iw==0 && myParams.bFix2Ddrift && nz==1) {
          timestamps_for_fitting = malloc(nphases * sizeof(float));
          refImage = malloc((nx+2)*ny*sizeof(float));
        }

      }

      for (direction=0; direction<ndirs; direction++) {

        /* data assumed taken with z changing fast, direction changing slowly */
        rawImages = (float **)  savedbands[direction];
        bands =  savedbands[direction];
        if (myParams.nOffImages) offImages = bufOffImages[direction];

        if (myParams.phaseSteps != 0) {
          // User specified the non-ideal (i.e., not 2pi/5) phase steps for each orientation
          // Under this circumstance, calculate the non-ideal sepMatrix:

          if (it==0 && iw==0 && direction==0)
            phaseList = malloc(nphases*sizeof(float));
          for (i=0; i<nphases; i++)
            phaseList[i] = i*myParams.phaseSteps[direction];

          makematrix(nphases, norders, direction, phaseList, sepMatrix, noiseVarFactors);
        }
        else if (myParams.phaseList)
          makematrix(nphases, norders, direction, myParams.phaseList+direction*nphases,
                     sepMatrix, noiseVarFactors);

        for (z=0; z<nz; z++) {
          int zsec;   /* to indicate which current section of the raw data file to read;
                         method varies depending on if it's fast SIM data or not. */
          /* First determine which section of the raw data to load */
          if (myParams.bFastSIM)  /* data organized into (nz, ndirs, nphases) */
            zsec = (z*nDirPhases + direction*nphases) * (1+myParams.nOffImages);
          else   /* data organized into (ndirs, nz, nphases) */
            zsec = (direction*nz*nphases + z*nphases) * (1+myParams.nOffImages);

          for (phase=0; phase<nphases; phase++) {

            /* If switch-off images are available, read them now */
            if (myParams.nOffImages) {
              load_and_flatfield(istream_no, zsec, iw, it, offImages[phase]+(z+zoffset)*nxy, buffer,
                                 nx, ny, background, backgroundExtra, slope, inscale, myParams.bUsecorr);
              for (i=1; i<myParams.nOffImages; i++) {
                load_and_flatfield(istream_no, zsec+i, iw, it, oneOffImage, buffer,
                                 nx, ny, background, backgroundExtra, slope, inscale, myParams.bUsecorr);
                image_arithmetic(offImages[phase]+(z+zoffset)*nxy, oneOffImage, nxy, 1., 1.);
              }
              zsec += myParams.nOffImages;
            }

            // do not average off images:
/*             if (myParams.nOffImages>1)  // average all the OFF images (may change later TODO) */
/*               image_arithmetic(offImages[phase]+z*nxy, NULL, nxy, 1./myParams.nOffImages, 0); */

            /* read the normal (or the only one in linear SIM) exposure*/
            if (myParams.bBgInExtHdr) { /* subtract the background value of each exposure stored in extended header */
              int extInts;
              float extFloats[3];
              IMRtExHdrZWT(istream_no, zsec, iw, it, &extInts, extFloats);
              /* for (i=0; i<nx*ny; i++) */
              /*   background[i] = extFloats[2]; */
              backgroundExtra = extFloats[2];
            }
            load_and_flatfield(istream_no, zsec, iw, it, rawImages[phase]+(z+zoffset)*nxy, buffer,
                               nx, ny, background, backgroundExtra, slope, inscale, myParams.bUsecorr);

            floatsection = rawImages[phase]+(z+zoffset)*nxy;

            if (myParams.napodize>=0) {
              apodize(myParams.napodize, nx, ny, floatsection);
              if (myParams.nOffImages)
                apodize(myParams.napodize, nx, ny, offImages[phase]+(z+zoffset)*nxy);
            }
            else if (myParams.napodize==-1) {
              cosapodize(nx, ny, floatsection);
              if (myParams.nOffImages)
                cosapodize(nx, ny, offImages[phase]+(z+zoffset)*nxy);
            }
            zsec ++;
          } /* end for (phase), loading, flatfielding, and apodizing raw images */


          if (myParams.do_rescale) {
            rescale(nx, ny, nz, z, zoffset, direction, iw, it, nphases, rawImages, myParams.equalizez, myParams.equalizet, sum_dir0_phase0);
            // also rescale offImages:
            if (myParams.nOffImages>1)
              rescale(nx, ny, nz, z, zoffset, direction, iw, it, nphases, offImages, myParams.equalizez, myParams.equalizet, sum_dir0_phase0_off);
          }

          /* 2D FFT all pre-processed, not yet separated, raw images */
          for (phase=0; phase<nphases; phase++) {
            fftwf_execute_dft_r2c(rfftplan, rawImages[phase]+(z+zoffset)*nxy, (fftwf_complex *) (rawImages[phase]+(z+zoffset)*nxy));
            if (myParams.nOffImages)
              fftwf_execute_dft_r2c(rfftplan, offImages[phase]+(z+zoffset)*nxy, (fftwf_complex *) (offImages[phase]+(z+zoffset)*nxy));
          }

          /* Estimate signal-to-noise ratio for each 2D section here, and use 1/SNR as the Wiener factor in filterbands() */
          /* For reference, see Wikipedia's Wiener deconvolution page */
          /* if (myParams.bUseEstimatedWiener) { */
          /*   myParams.wiener = sqrt(estimate_Wiener( (fftwf_complex**) rawImages, nx, ny, z+zoffset, nphases, rdistcutoff)); */
          /*   if (z == nz/2) printf("wiener=%f\n", myParams.wiener); */
          /* } */

          if (myParams.bFix2Ddrift) {
            /* standard separation matrix will be modified based on drift estimation */

            if (nz==1) { // otherwise, use the driftlist estimated above rigt after determinedrift_bt_dir()
              calcRefImage(rawImages, refImage,
                           offImages, myParams.nOffImages,
                           nx, ny, nphases,
                           myParams.type_of_refImage);

              determinedrift_2D((fftwf_complex**) rawImages, (fftwf_complex**) offImages, myParams.nOffImages, (fftwf_complex*)refImage,
                                drifts, nphases, nx, ny, direction, rdistcutoff, dkx, dky, myParams.drift_filter_fact);

              if (myParams.bFitDrift) {
                // obtain time stamps for fitting
                for (phase=0; phase<nphases; phase++)
                  timestamps_for_fitting[phase] = timestamps[direction*nphases*(1+myParams.nOffImages)
                                                             + phase*(1+myParams.nOffImages)
                                                             + myParams.nOffImages];
                fitXYdrift(drifts, timestamps_for_fitting, nphases, driftlist+direction*nphases, timestamps_for_fitting, nphases);
              }
              else
                for (phase=0; phase<nphases; phase++) {
                  driftlist[direction*nphases+phase] = drifts[phase];
/*                   driftlist[direction*nphases+phase].x = drifts[phase].x; */
/*                   driftlist[direction*nphases+phase].y = drifts[phase].y; */
                }

              printf("\nWithin-dir drift estimates and fit:\n");
              for (phase=0; phase<nphases; phase++)
                printf("%d: (%f, %f),   (%f, %f)\n", phase, drifts[phase].x, drifts[phase].y,
                       driftlist[direction*nphases+phase].x,  driftlist[direction*nphases+phase].y);
            }

            // for each exposure and switch-off image (if available), fix the 2D drift
            fixdrift_2D((fftwf_complex **) rawImages, driftlist, nphases, nx, ny, nz, direction, z);
            if (myParams.nOffImages)
              fixdrift_2D((fftwf_complex **) offImages, driftlist, nphases, nx, ny, nz, direction, z);

            /* Saved the aligned raw data if so requested */
            if (myParams.bSaveAlignedRaw) {
              float *temp = malloc(nxy*sizeof(float));
              fftwf_plan rfftplan2B = fftwf_plan_dft_c2r_2d(ny, nx, (fftwf_complex *) temp, temp, FFTW_ESTIMATE);
              for (phase=0; phase<nphases; phase++) {
                if (myParams.nOffImages) {
                  memcpy(temp, offImages[phase]+z*nxy, nxy*sizeof(float));
                  fftwf_execute(rfftplan2B);
                  for (i=0; i<ny; i++)
                    IMWrLin(aligned_stream_no, temp+i*(nx+2));
                }
                memcpy(temp, rawImages[phase]+z*nxy, nxy*sizeof(float));
                fftwf_execute(rfftplan2B);
                for (i=0; i<ny; i++)
                  IMWrLin(aligned_stream_no, temp+i*(nx+2));
              }
              free(temp);
              fftwf_destroy_plan(rfftplan2B);
            }

            if (myParams.bFixPhaseStep) {
              // makematrix() takes into account phase error, based on drift estimate and acquisition-time phase correction
              calcPhaseList(phaseList, driftlist, phaseAbs,
                            myParams.k0startangle+direction*delta_angle, myParams.linespacing, dxy,
                            myParams.nOffImages, nphases, nz, direction, z);

              makematrix(nphases, norders, direction, phaseList, sepMatrix, noiseVarFactors);
            }
          }

          separate(nx, ny, z+zoffset, direction, nphases, norders, rawImages, sepMatrix);  // unmixing of components in real or reciprocal space

          if (myParams.nOffImages) {
            separate(nx, ny, z+zoffset, direction, nphases, norders, offImages, sepMatrix);
            if (myParams.bCombineOffExp)
              combine_exp_off_F(bands, (fftwf_complex **) offImages, nx, ny, norders, direction, rdistcutoff,
                                dkx, dky, noiseVarFactors, myParams.combineOffOrders);
          }


        }    /* end for (z) */  /* reading in actual data */


        if (nz>1)  // 1D FFT of a stack of 2D FFTs
          for (i=0; i<nphases; i++)
            fftwf_execute_dft(fftplan1D, bands[i], bands[i]);

        /* maybe estimate signal-to-noise here after 3D FFT is done? */
        if (myParams.bUseEstimatedWiener) {
          est_wiener[direction] = estimate_Wiener_3D(bands[0], nx, ny, nz, dkx, dky, dz, wave[0], &myParams);
        }

        if (myParams.bMakemodel)
          /* Use the OTF to simulate an ideal point source; replace bands with simulated data */
          makemodeldata(nx, ny, nz0, bands, norders, k0[direction], dkx, dky, dz, otf, wave[iw], &myParams);

        /* save intermediate results if requested */
        if (myParams.bSaveSeparated) {
          for (phase=0; phase<nphases; phase++)
            for (z=0; z<nz; z++)
              IMWrSec(separated_stream_no, rawImages[phase]+(z+zoffset)*nxy);
          continue; // skip the k0 search and modamp fitting
        }

        /* After separation and FFT, the array rawImages, now referred to as the array bands[], contains the center
           band (bands[0]), and the real (bands[1], bands[3], etc.) and imaginary part (bands[2], bands[4], etc.) bands
           of the different orders */

        if (! (ntimes>1 &&  it>0 && myParams.bUseTime0k0)) {
          k0[direction] = k0guess[direction];
          printf("k0guess[direction %d] = (%f, %f) 1/micron\n", direction, k0guess[direction].x, k0guess[direction].y);
        }

        /* Now to fix 3D drift between dirs estimated by determinedrift_3D() */
        if (direction!=0 && myParams.bFixdrift)
          fixdrift_bt_dirs(bands, norders, drift_bt_dirs[direction], nx, ny, nz0);

        if (myParams.bSearchforvector &&   /* assume k0 vector not well known, so fit for it */
            !(ntimes>1 && it>0 && myParams.bUseTime0k0)) { /* In time series, can choose to use the time 0 k0 fit for the rest of the series */
          /* Find initial estimate of modulation wave vector k0 by cross-correlation. */
          findk0(bands, overlap0, overlap1, nx, ny, nz0, norders, &(k0[direction]), dxy, dz, otf, wave[iw], &myParams);

          if (myParams.bSaveOverlaps) // output the overlaps
            for (z=0; z<nz0; z++) {
              IMWrSec(overlaps_stream_no, overlap0+z*nx*ny);
              IMWrSec(overlaps_stream_no, overlap1+z*nx*ny);
            }

          printf("Initial guess by findk0() of k0[direction %d] = (%f,%f) 1/micron\n", direction, k0[direction].x, k0[direction].y);

          /* refine the cross-corr estimate of k0 by a search for best k0 vector direction
             and magnitude using real space waves*/

          printf("before fitk0andmodamp\n");
          /* if (myParams.recalcarrays==0 && dist>1.0) */ myParams.recalcarrays = 1;
          /* if k0 is very close to the guess, we can save time by not recalculating the overlap arrays */

          fitk0andmodamps(bands, overlap0, overlap1, nx, ny, nz0, norders, &(k0[direction]),
                          dxy, dz, otf, wave[iw], amp[direction], &myParams);

          if (it==0)
            k0_time0[direction] = k0[direction];
          /* check if the k0 vector found is reasonably close to the guess */
          deltak0.x = k0[direction].x - k0guess[direction].x;
          deltak0.y = k0[direction].y - k0guess[direction].y;
          dist = sqrt(deltak0.x * deltak0.x + deltak0.y * deltak0.y);
          if (dist/k0magguess > K0_WARNING_THRESH)  printf("WARNING: ");
          printf("best fit for k0 is %f percent off expected value.\n", dist/k0magguess*100);
          if (ntimes>1 && it>0 && myParams.bUseTime0k0 && dist > 2*K0_WARNING_THRESH) {
            k0[direction] = k0_time0[direction];
            printf("k0 estimate of time point 0 is used instead\n");
            for (order=1; order<norders; order++) {
              if (nz0>1)
                corr_coeff = findrealspacemodamp(bands, overlap0, overlap1, nx, ny, nz0, 0, order, k0[direction], dxy, dz,
                                                 otf, wave[iw], &amp[direction][order], &amp_inv, &amp_combo, 1, &myParams);
              else
                corr_coeff = findrealspacemodamp(bands, overlap0, overlap1, nx, ny, nz0, order-1, order, k0[direction], dxy, dz,
                                                 otf, wave[iw], &amp[direction][order], &amp_inv, &amp_combo, 1, &myParams);
              printf("modamp mag=%f, phase=%f\n",cabsf(amp[direction][order]), cargf(amp[direction][order]));
            }
          }

        }
        else {  /* assume k0 vector known, so just fit for the modulation amplitude and phase */

          printf("known k0 for direction %d = (%f, %f) \n", direction, k0[direction].x, k0[direction].y);
          for (order=1; order<norders; order++) {
            corr_coeff = findrealspacemodamp(bands, overlap0, overlap1, nx, ny, nz0, 0, order, k0[direction], dxy, dz,
                                             otf, wave[iw], &amp[direction][order], &amp_inv, &amp_combo, 1, &myParams);
            printf("modamp mag=%f, phase=%f\n",cabsf(amp[direction][order]), cargf(amp[direction][order]));
            printf("reverse modamp mag=%f, phase=%f\n", 1/cabsf(amp_inv), -cargf(amp_inv));
            printf("combined modamp mag=%f, phase=%f\n", cabsf(amp_combo), cargf(amp_combo));
			printf("correlation coeff=%f\n\n", corr_coeff);
            if (order==1 && myParams.bSaveOverlaps) // output the overlaps
              for (z=0; z<nz0; z++) {
                IMWrSec(overlaps_stream_no, overlap0+z*nx*ny);
                IMWrSec(overlaps_stream_no, overlap1+z*nx*ny);
              }

          }
        }     /* if(searchforvector) ... else ... */

        if (nz==1) /* In 2D SIM, amp stores modamp's between each adjacent pair of bands. We want to convert this to modamp w.r.t. order 0 */
          for (order=2; order < norders; order++)
            amp[direction][order] *= amp[direction][order-1];

        if (myParams.forceamp[1] > 0.0)   /* force modamp's amplitude to be a value user provided (ideally should be 1)  */
          for (order=1;order<norders;order++) {
            float a = cabsf(amp[direction][order]);
            if (a<myParams.forceamp[order]) {
              ampfact = myParams.forceamp[order]/a;
              amp[direction][order] *= ampfact;
              printf("modamp mag=%f, phase=%f  \n",cabsf(amp[direction][order]),
                     cargf(amp[direction][order]));
            }
          }

        /* In 2D NLSIM, we often don't trust the modamp fit between neighboring high-order components; */
        /* only if fitallphases is True do all fitted modamps get actually used; */
        /* otherwise, order 2 and above's global phase is inferred from order 1's phase.  */
        if (!myParams.bFitallphases) {
          float base_phase = cargf(amp[direction][1]);
          fftwf_complex expiphi,  amplitude;
          float phi;
          for(order=2; order<norders; order++) {
            amplitude = cabsf(amp[direction][order]) + 0*I;
            phi = order * base_phase /*+ M_PI*(order-1)*/ ;   /* sign flip for every other order happens only in saturation case? */
            expiphi = cos(phi) + I * sin(phi);
            amp[direction][order] = amplitude * expiphi;
          }
        }
        printf("\nk0[%i].x=%f, k0[%i].y=%f\n", direction, k0[direction].x, direction, k0[direction].y);

      } /* end for (direction)  */   /* done finding the modulation vectors and phases for all directions */

        
      if (myParams.bUseEstimatedWiener) {
        myParams.wiener = sum_of_pixels(est_wiener, ndirs, 1, 0)/ndirs;
        printf("\nEstimated wiener average is %f\n\n", myParams.wiener);
      }

      if (it == ntimes-1) {
        if (myParams.bSaveSeparated) {
          IMWrHdr(separated_stream_no, "separated bands of all directions", 1, 0, 1, 0);
          IMClose(separated_stream_no);
        }
        if (myParams.bSaveAlignedRaw) {
          IMWrHdr(aligned_stream_no, "drift-corrected raw images", 1, aligned_header.amin, aligned_header.amax, aligned_header.amean);
          IMClose(aligned_stream_no);
        }
        if (myParams.bSaveOverlaps) {
          IMWrHdr(overlaps_stream_no, "overlaps in real space", 1, 0, 1, 0);
          IMClose(overlaps_stream_no);
        }
      }

      if (myParams.bSaveSeparated || myParams.bSaveAlignedRaw || myParams.bSaveOverlaps) {
        if (it < ntimes-1)
          continue;
        else {
          printf("\nQuit processing because either bSaveSeparated, bSaveAlignedRaw, or bSaveOverlaps is TRUE\n");
          return 0;
        }
      }

      /* Now that we know the modulation vectors and amplitudes, start assembling the final data set, one direction at a time. */

      if (it == 0 && iw == 0) {
        /* bigbuffer is work place for preparing for the bands assembly */
        bigbuffer  = (fftwf_complex *) fftwf_malloc((unsigned) ((myParams.zoomfact*nx) *
                                                                 (myParams.zoomfact*ny) *
                                                                 (myParams.z_zoom*nz) * sizeof(fftwf_complex)
                                                                )
                                                    );
        /* outbuffer holds the final assembly */
        outbuffer  = (float *) malloc((unsigned)((myParams.zoomfact*nx)*(myParams.zoomfact*ny)*(myParams.z_zoom*nz0)*sizeof(float)));
        if (!bigbuffer || !outbuffer) {
          printf("Not enough memory for bigbuffer or outbuffer\n");
          return -1;
        }
        if (nz0 > 1)
          fftplan  = fftwf_plan_dft_3d(myParams.z_zoom*nz0, (int) (myParams.zoomfact*ny), (int) (myParams.zoomfact*nx), bigbuffer, bigbuffer, FFTW_BACKWARD, FFTW_ESTIMATE);
        else
          fftplan = fftwf_plan_dft_2d((int) (myParams.zoomfact*ny), (int) (myParams.zoomfact*nx), bigbuffer, bigbuffer, FFTW_BACKWARD, FFTW_ESTIMATE);
      }

      memset(outbuffer, 0, (unsigned)((myParams.zoomfact*nx)*(myParams.zoomfact*ny)*(myParams.z_zoom*nz0)*sizeof(float)));

      if (myParams.bPreciseApoBoundary) {
        // calculate angular boundaries between sectors defined by the highest-harmonics' arcs (the "petals")
        boundaryAngles = malloc(ndirs * 2 * sizeof(float));
        calculateSectorBoundaries(ndirs, k0, boundaryAngles);
      }

      for (direction=0; direction<ndirs; direction++) {
        printf("for direction %d, do filterbands() and assemblerealspace()\n", direction);

        bands=savedbands[direction];
        printf("before filterbands\n");

        /* Prepare the bands to be assembled: Apply appropriate OTF-based weighting in the overlaps,
           incorporating a wiener-like filter  */
        filterbands(direction, bands, k0, ndirs, norders, otf, dxy, dz,
                    amp, noiseVarFactors, nx, ny, nz0, wave[iw], &myParams);

        /* FFT back to real space and assemble */
        printf("before assemblerealspacebands\n");
        assemblerealspacebands(direction, outbuffer, bigbuffer, bands, ndirs, norders, k0, dxy, nx, ny, nz0,
                               fftplan, myParams.zoomfact, myParams.z_zoom, myParams.explodefact, myParams.bNoOrder0);

      } /* end for (dir) */

      /* Write out  */

      printf("before write\n");
      {
        float *ptr = outbuffer + (int)(zoffset*myParams.z_zoom * nx*ny*myParams.zoomfact*myParams.zoomfact);
        for (i=0; i<nz*myParams.z_zoom; i++) {
          IMWrSec(ostream_no, ptr);
          if (it==0) {
            for (j=0; j<(int) (nx*ny*myParams.zoomfact*myParams.zoomfact); j++)
              if (ptr[j] > maxval)
                maxval = ptr[j];
              else if (ptr[j] < minval)
                minval = ptr[j];
          }
          ptr += (int) (myParams.zoomfact*nx*myParams.zoomfact*ny);
        }
        if (it==0 && iw==0) {
          header.amin = minval;
          header.amax = maxval;
        }
      }

      printf("Time point %d, wave %d done\n",it, iw);

      /* if (ntimes==1 && nwaves==1) */
        /* fftwf_destroy_plan(fftplan); /\* Don't destroy it if it will be reused *\/ */

    } /* end for (wave) */
  } /* end for (it)  (loop over time points) */
  
  IMClose(istream_no);
  
  fftwf_destroy_plan(rfftplan);
  fftwf_destroy_plan(fftplan);
  fftwf_free(bigbuffer);
  free(outbuffer);

  /* Save the commandline into the header */
  titles[0] = '\0';
  for (i=3; i<argc; i++) {
    strcat(titles, argv[i]);
    strcat(titles, " ");
  }
  IMAlLab(ostream_no, titles, strlen(titles)/80+1);
  IMWrHdr(ostream_no, header.label, 1, header.amin, header.amax, header.amean);
  IMClose(ostream_no);

  free(sum_dir0_phase0);
  free(sum_dir0_phase0_off);
  free(k0);
  free(k0_time0);
  free(k0guess);
  free(drift_bt_dirs);
  free(sepMatrix);
  free(noiseVarFactors);
  for (direction=0; direction<ndirs; direction++) {
    free(amp[direction]);
    for (i=0; i<nphases; i++) {
      /* free(savedbands[direction][i]); */
      if (myParams.nOffImages)
        free(bufOffImages[direction][i]);
    }
    free(savedbands[direction]);
    free(bufOffImages[direction]);
  }
  free(amp);
  free(savedbands);
  free(bufOffImages);

  fftwf_free(overlap0);
  fftwf_free(overlap1);

  free(buffer);
  for (i=0; i<norders; i++)
    free(otf[i]);
  free(otf);
  free(background);
  free(slope);
  if (expDose) free(expDose);
  if (timestamps) free(timestamps);
  if (phaseAbs) free(phaseAbs);

  return 0;
}


/*********************************   apodize   **********************************************/
/*   softens the edges of a singe xy section to reduce edge artifacts and improve the fits  */
/********************************************************************************************/
void apodize(int napodize,int nx,int ny,float *image)
{
  float diff,fact;
  int k,l;

  for (k=0;k<nx;k++) {
    diff = (image[ (ny-1)*(nx+2) +k ] - image[ /* 0*(nx+2)+ */ k ]) / 2;
    for (l=0;l<napodize;l++) {
      fact = 1 - sin((((float)l+0.5)/napodize)*M_PI*0.5);
      image[l*(nx+2) + k] += diff*fact;
      image[(ny-1-l)*(nx+2) + k] -= diff*fact;
    }
  }
  for (l=0;l<ny;l++) {
    diff = (image[ l*(nx+2) + nx-1 ] - image[ l*(nx+2) /* +0 */ ]) / 2;
    for (k=0;k<napodize;k++) {
      fact = 1 - sin((((float)k+0.5)/napodize)*M_PI*0.5);
      image[l*(nx+2) + k] += diff*fact;
      image[l*(nx+2) + nx-1-k] -= diff*fact;
    }
  }
}


/*******************************   cosapodize   ********************************/
/*   softens the edges to reduce edge artifacts and improve the fits           */
/*******************************************************************************/
void cosapodize(int nx,int ny,float *image)
{
  float xfact,yfact;
  int k,l;

  /*printf("in cosapodize\n");*/
  for (k=0;k<nx;k++) {
    xfact = sin(M_PI*((float)k+0.5)/nx);
    for (l=0;l<ny;l++) {
      yfact = sin(M_PI*((float)l+0.5)/ny);
      image[l*(nx+2) + k] *= (xfact*yfact);
    }
  }
}

/***************************** rescale *********************************/
/* Corrects for imperfectly equal exposure in the input images.        */
/* Does this by equalizing the sum.  This assumes that the             */
/* strength of the sample structure spatial frequency component that   */
/* exactly equals the modulation frequency is small compared to the    */
/* total intensity.  In some test samples, such as single beads,       */
/* that assumption may not be true.                                    */
/* When equalizez is 0, equalizes all the average values (all phases   */
/* and directions) for each z section, but does not equalize           */
/* between different z sections. Otherwise, equalize across z as well. */
/***********************************************************************/
void rescale(int nx, int ny, int nz, int z, int zoffset, int direction, int wave, int t, int nphases, float *rawImages[], int equalizez, int equalizet, double *sum_dir0_phase0)
/* zoffset -- used when there's z padding and rawImages is allocated more z sections than nz */
{
  double *sum;
  float ref;
  int k, l, phase, nxy, rowstart;

  nxy=(nx+2)*ny;
  sum = (double *)malloc(nphases * sizeof(double));
#pragma omp parallel for private(phase, l, k, rowstart)
  for (phase=0;phase<nphases;phase++) {
    sum[phase] = 0.0;
    for (l=0;l<ny;l++) {
      rowstart = (z+zoffset)*nxy + l*(nx+2);
      for (k=0;k<nx;k++)
        sum[phase] += rawImages[phase][rowstart + k];
    }
  }
  if (direction==0 && !(equalizet && t!=0)) /* first direction for this z */  /* if equalizet is true, equalize all time points in the same way as the first one */
    sum_dir0_phase0[wave*nz + z] = sum[0];           /* for later directions, normalize to the first direction */

  if (equalizez)   /* equalize all z sections as well as the phases and directions within each z section */
    ref = sum_dir0_phase0[wave*nz+0];
  else           /* only equalize the phases and directions within each z section */
    ref = sum_dir0_phase0[wave*nz+z];

#pragma omp parallel for private(phase, l, k, rowstart)
  for (phase=0;phase<nphases;phase++) {
    float ratio;
    ratio= ref/sum[phase];
    for (l=0;l<ny;l++) {
      rowstart = (z+zoffset)*nxy + l*(nx+2);
      for (k=0;k<nx;k++)
        rawImages[phase][rowstart + k] *= ratio;
    }
  }

  free(sum);
}


/***************************** makematrix ************************************/
/*     generates the matrix that is to be used to separate the raw indata    */
/*     into the different bands of sample information.                       */
/*  Two cases:
    1. If arrPhases is NULL, we're certain that phases are equally spaced within 2Pi range;
    2. Else, user supplies a list of phases actually used (e.g, when there is drift or 
    non-ideal SLM patterns are used), and noise variance amplification factors for each 
    order and direction is calculated based on the non-ideal separation matrix
*/
/*****************************************************************************/
void makematrix (int nphases, int norders, int dir, float *arrPhases, float *sepMatrix, float * noiseVarFactors)
/*
  norders is not necessarily nphases/2+1, in cases such as nonlinear SIM where the number of output orders can be smaller.
 */
{
  int j, order;
  float phi;

  // norders could be less than (nphases+1)/2
  if (arrPhases == 0) { /* phase values equally spaced between 0 and 2Pi */
    phi = 2*M_PI/nphases;
    for (j=0; j<nphases; j++) {
      sepMatrix[0*nphases+j] = 1.0;
      for (order=1; order<norders; order++) {
        sepMatrix[(2*order-1)*nphases+j] = cos(j*order*phi);  /* With this definition, bandplus = bandre + i bandim has coefficient exp() with unit normalization */
        sepMatrix[2*order    *nphases+j] = sin(j*order*phi);
      }
    }
  }
  else {/* User supplied phase values in arrPhases */
    int *ipvt, info, i, nCols;
    float *workspace, *forwardM, *A_t_A, unit=1., zero=0.;

    /* Here we differentiate between 2 options: */
    /* 1. use direct matrix inversion;          */
    /* 2. use least-square solution.            */

    nCols = 2*norders-1;

    if (nCols < nphases) {  // use least-square
      ipvt = malloc(nCols*sizeof(int));
      forwardM   = malloc(nphases*nCols*sizeof(float));
      A_t_A = malloc(nCols*nCols*sizeof(float));
      /* First construct the forward matrix*/
      for (i=0; i<nphases; i++) {
        forwardM[i] = 1.0 / nphases;
        for (order=1; order<norders; order++) {
          forwardM[i+(2*order-1)*nphases] = 2 * cos(order*arrPhases[i]) / nphases;
          forwardM[i+2*order*nphases]     = 2 * sin(order*arrPhases[i]) / nphases;
        }
      }

/*       print_matrix("forward matrix:", nphases, nCols, forwardM, nphases); */

      /* multiply transposed forward matrix with forward matrix */
      sgemm_("T", "N", &nCols, &nCols, &nphases, &unit, forwardM,  &nphases, forwardM, &nphases,  &zero, A_t_A,  &nCols);
/*       print_matrix("A transpose times A:", nCols, nCols, A_t_A, nCols); */

      /* Then invert the forward matrix to form the sep matrix*/
      sgetrf_(&nCols, &nCols, A_t_A, &nCols, ipvt, &info);
      if (info==0) { // successful factorization
        workspace = malloc(nCols*sizeof(float));
        sgetri_(&nCols, A_t_A, &nCols, ipvt, workspace, &nCols, &info);
        free(workspace);
      }
      else
        printf("WARNING: sgetri() returns non 0.\n");

      if (info!=0)
        printf("WARNING: sgetrf() returns non 0.\n");

      // multiply inverse A_t_A with transposed forward matrix
      sgemm_("N", "T", &nCols, &nphases,  &nCols, &unit, A_t_A,  &nCols, forwardM, &nphases,  &zero, sepMatrix, &nCols);

/*       print_matrix("Separation Matrix:", nCols, nphases, sepMatrix, nCols); */
      /* transpose back to C-style matrix */
      matrix_transpose(sepMatrix, nphases, nCols);

      free(ipvt);
      free(forwardM);
      free(A_t_A);
    }
    else {  // use direct inversion
      ipvt = malloc(nphases*sizeof(int));
      /* First construct the forward matrix (in C-style storage convention ) */
      for (i=0; i<nphases; i++) {
        sepMatrix[i*nphases] = 1.0 / nphases;
        for (order=1; order<norders; order++) {
          sepMatrix[i*nphases+(2*order-1)] = 2 * cos(order*arrPhases[i]) / nphases;
          sepMatrix[i*nphases+2*order]     = 2 * sin(order*arrPhases[i]) / nphases;
        }
      }
      /* Then invert the forward matrix to form the sep matrix*/
      sgetrf_(&nphases, &nphases, sepMatrix, &nphases, ipvt, &info);
      if (info==0) { // successful factorization
        workspace = malloc(nphases*sizeof(float));
        sgetri_(&nphases, sepMatrix, &nphases, ipvt, workspace, &nphases, &info);
        free(workspace);
      }
      else
        printf("WARNING: sgetri() returns non 0.\n");

      if (info!=0)
        printf("WARNING: sgetrf() returns non 0.\n");
      free(ipvt);
    }

    /* Report noise factors */
    for (order=0;order<norders;order++)   noiseVarFactors[dir*norders+order] = 0.0;
    for (j=0;j<nphases;j++) {
      noiseVarFactors[dir*norders+0] += pow(sepMatrix[0*nphases+j], 2);
      for (order=1;order<norders;order++)
        noiseVarFactors[dir*norders+order] += pow(sepMatrix[(2*order-1)*nphases+j], 2) + pow(sepMatrix[2*order*nphases+j], 2);
    }

    printf(" Separation noise factors: ");
    for (order=0;order<norders;order++) {
      noiseVarFactors[dir*norders+order] /= nphases;
      printf(" Order %d: %.3f,",order, sqrt(noiseVarFactors[dir*norders+order]));
    }
    printf("\n");

  }

  printf("Separation matrix:\n");
  for (j=0;j<norders*2-1;j++) {
    int i;
    for (i=0; i<nphases; i++)
      printf("%9.5f ", sepMatrix[j*nphases+i]);
    printf("\n");
  }
  printf("\n");
}


/*****************************  findk0  ***********************************/
/*     Finds the modulation wave vector by a cross-correlation of the     */
/*     different bands in the region where they overlap (and where they   */
/*     should in principle only differ by a constant complex factor).     */
/**************************************************************************/

void findk0(fftwf_complex *bands[], fftwf_complex *overlap0, fftwf_complex *overlap1, int nx, int ny, int nz, int norders, vector *k0, float dxy, float dz, fftwf_complex *OTF[], short wave, ReconParams * pParams)
{
  fftwf_complex Xval, Yval;
  int i, j, z, ind, xyind, nxy;
  int  fitorder1,  fitorder2;
  float *crosscorr;
  fftwf_complex prod, *crosscorr_c;
  fftwf_plan fftplan2d;
  vector old_k0;
  float dkx, dky;

  crosscorr  =  (float *) malloc((nx*ny) * sizeof(float));
  crosscorr_c = (fftwf_complex *) calloc((nx*ny), sizeof(fftwf_complex));
  fftplan2d = fftwf_plan_dft_2d(ny, nx, crosscorr_c, crosscorr_c, FFTW_FORWARD, FFTW_ESTIMATE);

  /* make arrays that contain only the overlapping parts of fourier space */
  /* otf-compensate in that region, set to zero elsewhere. Transform to real space.  */

  /* For mats's non-linear saturation experiment, we can't use order 0 and norders-1 to fit k0 since there's no overlap between them */
  /* The idea is to change makeoverlaps() so that it can make overlaps between any pair of bands, intra- or inter-directions */
  /* For now, the inter-direction overlap is not implemented because of the way we handle the input -- read in one direction of raw data
     at a time so as not to overwhelm memory */
  fitorder1 = 0;
  if (nz>1) {
    if (!pParams->bBessel)
      fitorder2 = 2;
    else
      fitorder2 = 1; //norders-1;
  }
  else
    fitorder2 = 1;

  // k0 is in 1/micron
  makeoverlaps(bands, overlap0, overlap1, nx, ny, nz, fitorder1, fitorder2, (*k0).x, (*k0).y, dxy, dz, OTF, wave, pParams);

  /* the returned overlap arrays are in real space, but are still complex since we cut out a non-hermitian part in fourier space */

  /*  multiply and retransform (i.e. cross-correlate).  Average over z to get 2D lateral-only correlation *****/
  for (z=0;z<nz;z++)
#pragma omp parallel for private(i,j,ind,xyind,Xval,Yval,prod) shared(crosscorr_c)
    for (i=0;i<ny;i++)
      for (j=0;j<nx;j++) {
        ind = z*nx*ny+i*nx+j;
        xyind = i*nx+j;
        Xval = overlap0[ind];
        Yval = overlap1[ind];
        prod = Xval * conjf(Yval);
        crosscorr_c[ xyind ] += prod;
      }
  fftwf_execute(fftplan2d);

  nxy = nx*ny;
#pragma omp parallel for private(i) shared(crosscorr)
  for (i=0;i<nxy;i++)
    crosscorr[i] = mag2(crosscorr_c[i]);
  /******** end cross-correlation *************/

  old_k0 = *k0;

  findpeak(crosscorr, nx, ny, k0);  /* position of the cross-correlation peak is returned in the vector k0  */
  // k0 returned in by findpeak() is in pixels; needs to be converted back to 1/micron

  /* Lin 11/25/02: Let's try the following instead for the time being. (introduced old_k0 for this purpose) */
  /* if((*k0).y > ny/2 && old_k0.y < 0) (*k0).y -= ny;
     if((*k0).x > nx/2 && old_k0.x < 0) (*k0).x -= nx; */
  /* Mats 01/16/2008: I think Lin's concern above is that the k0 position can appear translated by nx or ny
     due to aliasing if the camera undersamples relative to the pattern frequency. I think all one has to do is to
     make sure that the new k0 is closer to the old (guessed) k0 than to k0+nx or k0-nx. */

  dkx = 1 / (nx*dxy);
  dky = 1 / (ny*dxy);
  if (old_k0.x * dkx < (*k0).x-nx/2) (*k0).x -= nx;
  if (old_k0.x * dkx > (*k0).x+nx/2) (*k0).x += nx;
  if (old_k0.y * dky < (*k0).y-ny/2) (*k0).y -= ny;
  if (old_k0.y * dky > (*k0).y+ny/2) (*k0).y += ny;

  k0->x *= dkx/fitorder2;
  k0->y *= dky/fitorder2; /* return k0 (in 1/micron) of the first order regardless of fitorder2 */

  free(crosscorr);
  free(crosscorr_c);
  fftwf_destroy_plan(fftplan2d);
  return;
}


/*********************************  makeoverlaps  *********************************/
/*    Makes arrays ...part that contain only those part of the input bands        */
/*    which in fourier space overlaps the other band which is offset by (k0x,k0y).*/
/*    Filters the overlap region such that the information there ends up having   */
/*    an effective otf value (scale factor) = OTF1*OTF2/sqrt(OTF1^2 + OTF2^2).    */
/*    In another call to this function, the other band will get the same scaling, */
/*    so that the two outarrays should in theory differ only by a constant        */
/*    complex factor (the modulation amplitude).                                  */
/**********************************************************************************/

void makeoverlaps(fftwf_complex *bands[], fftwf_complex *overlap0, fftwf_complex *overlap1, int nx, int ny, int nz, int order1, int order2, float k0x, float k0y, float dxy, float dz, fftwf_complex *OTF[], short wave, ReconParams * pParams)
{
  int nxy, nxyin, nxyz, i, j, iin, jin, z0, z, xcent, ycent, indin, indout, conj;
  float rdistcutoff, zdistcutoff/* , zdistcutoff_order0, zdistcutoff_orderhighest */;
  fftwf_complex val1re, val1im, val2re, val2im, *band1re, *band1im, *band2re, *band2im;
  float root;
  float dkx, dky, dkz, kzscale, kx, ky, otfcutoff=0.006 /*.004*/;
  float rdist1, rdist12, rdist21, x1, y1, x12, y12, x21, y21;
  float lambdaem, alpha, k0mag;
  fftwf_complex otf1, otf2, otf21, otf12, fact;
  fftwf_plan fftplan;
  /* When using uniformly scaled OTF, the factor between order 0 and 2's magnitude.
     we need it because otfcutoff is used to decide whether to include a point in overlaps */
  float order0_2_factor=1.0;  /* only has other value if it's 3D case because of the OTF magnitude difference */


  if (nz>1) {
    if (!pParams->bBessel)
      order0_2_factor = 5;
    else
      order0_2_factor = 20;
  }

  xcent = nx/2;  ycent = ny/2;
  nxy = nx*ny; nxyin = (nx/2+1)*ny;
  nxyz = nx*ny*nz;
  dkx = 1/(nx*dxy);  /* inverse microns per pixel in data */
  dky = 1/(ny*dxy);
  if (dz>0)     /* sometimes in 2D SI data, dz is set to 0, which will mess up kzscale and otfinterpolate() later */
    dkz = (1/(nz*dz));   /* inverse microns per pixel in data */
  else
    dkz = pParams->dkzotf;
  kzscale = dkz / pParams->dkzotf;   /* ratio of axial direction pixel scales of data and otf */
  rdistcutoff = pParams->na*2/(wave/1000.0);    /* OTF support radial limit in 1/micron */
  if (rdistcutoff > 1/(2*dxy)) rdistcutoff = 1/(2*dxy);

  /* k0pix =  sqrt(k0x*k0x + k0y*k0y);   /\* k0 magnitude (for highest order) in pixels *\/ */
  k0mag = sqrt(k0x*k0x + k0y*k0y);   /* k0 magnitude (for lowest order) in 1/micron */
  lambdaem = (wave/pParams->nimm)/1000.0;  /* emission wavelength in the sample, in microns */
  alpha = asin(pParams->na/pParams->nimm);  /* aperture angle of objectives */

  /* zdistcutoff depends on options -- single lens or double lenses or light-sheet */
  if (!pParams->bTwolens & !pParams->bBessel)
    zdistcutoff = (int) ceil(((1-cos(alpha))/lambdaem) / dkz);    /* OTF support axial limit in data pixels */
  else if (pParams->bBessel) {  /* use the smaller kz cutoff of the two overlapping orders */
    float kzExMax, halfangle;

    kzExMax = 2*BesselNA / BesselLambdaEx;
    /* if (order1 == 2 || order2 == 2) { */
    /*   halfangle = acos(k0mag*2 / kzExMax); */
    /*   zdistcutoff =ceil( (kzExMax * sin(halfangle) + (1-cos(alpha))/lambdaem) / dkz); */
    /* } */
    /* else if (order1==1 || order2 == 1) { */
    /*   halfangle = acos(k0mag / kzExMax); */
    /*   zdistcutoff = ceil((kzExMax * sin(halfangle) + (1-cos(alpha))/lambdaem) / dkz); */
    /* } */
   halfangle = acos(k0mag* order2 / kzExMax); 
   zdistcutoff =ceil( (kzExMax * sin(halfangle) + (1-cos(alpha))/lambdaem) / dkz);
  }
  else {
    fprintf(stderr, "Sorry, this program doesn't handel 2-objective mode data\n");
    exit(-1);
  }

  if (zdistcutoff>nz/2)  zdistcutoff=((nz/2-1) > 0 ? (nz/2-1) : 0); // so that in case of 2D data, zdistcutoff is 0


  if (order1 == 0) {
    band1re = bands[0];
  }
  else {
    band1re = bands[order1*2-1];
    band1im = bands[order1*2];
  }

  /* Shall we assume that order 0 should never appear in the second fit order, i.e., order2 is never 0? */
  band2re = bands[order2*2-1];
  band2im = bands[order2*2];

  /* make arrays that contain only the overlapping parts of fourier space */
  /* otf-compensate in that region, set to zero elsewhere  */

  if (nz>1)
    fftplan = fftwf_plan_dft_3d(nz, ny, nx, overlap0, overlap0, FFTW_BACKWARD, FFTW_ESTIMATE);
  else
    fftplan = fftwf_plan_dft_2d(ny, nx, overlap0, overlap0, FFTW_BACKWARD, FFTW_ESTIMATE);

#pragma omp parallel for private(i) shared(overlap0, overlap1)
  for(i=0;i<nxyz;i++) {   /* have to re-zero entire array since it has been transformed to real space */
    overlap0[i] = 0;
    overlap1[i] = 0;
  }

  kx = k0x*(order2-order1);
  ky = k0y*(order2-order1);

#pragma omp parallel for private(i, j, x1, y1, conj, otf1,otf2,otf21,otf12, root, fact, val1re, val1im, val2re, val2im, iin, jin, z0, z, indin, indout, rdist1, rdist12, rdist21, x12, y12, x21, y21) shared(overlap0, overlap1)
  for (i=0;i<ny;i++)  /* coords of the output array. Full complex fft array, with f.sp. origin at i=j=0 */
    for (j=0;j<nx;j++) {
      x1 = j; y1 = i;
      if (x1>xcent) x1-=nx; /* fourier space coords with origin at x1=y1=0 and going  betwen -xcent and xcent */
      if (y1>ycent) y1-=ny;

      // Convert pixels into 1/micron:
      x1 *= dkx;
      y1 *= dky;

      rdist1 = sqrt(x1*x1+y1*y1);   /* distance from origin of input data */
      if (rdist1<=rdistcutoff) {     /* is x1,y1 within the theoretical lateral OTF support of
                                         the data that is to be scaled? */
        x12 = x1-kx; y12 = y1-ky;    /* f. sp. coords of the offset side band frame, if x1,y1 refer to the absolute frame */
        /* (kx,ky) is the center of the shifted sideband viewed from the absolute frame */
        rdist12 = sqrt(x12*x12+y12*y12);
        x21 = x1+kx; y21 = y1+ky;    /* the absolute f. sp. coords, if x1,y1 refer to the offset side band frame */
        /* (-kx,-ky) is the center of the center band viewed from the sideband frame */
        rdist21 = sqrt(x21*x21+y21*y21);
        if (rdist12<=rdistcutoff || rdist21<=rdistcutoff) {  /* within either theoretical OTF overlap, laterally */
          if (j<=xcent) {
            iin = i;   /* coords of input arrays */
            jin = j;
            conj = 0;
          }
          else {
            jin = nx-j;
            iin = (ny-i)%ny;     // "%ny" is to deal with i=0 !!!
            conj = 1;
          }
          if (rdist12<=rdistcutoff) { // within theoretical OTF support of the shifted side band

            for (z0=-zdistcutoff; z0<=zdistcutoff; z0++) {
              if (!(z0==0 && pParams->bNoKz0)) { /* use kz==0 plane or not? */
                otf1 = otfinterpolate(OTF[order1], x1, y1, z0, kzscale, pParams);
                if (cabsf(otf1)>otfcutoff) {     /* within the area of usably large OTF values of the data that is to be scaled */

                  otf12 = otfinterpolate(OTF[order2], x12, y12, z0, kzscale, pParams);
                  if (cabsf(otf12)*order0_2_factor>otfcutoff) {  /* within usable OTF support of the shifted side band */
                    if (conj) z = -z0;
                    else z = z0;
                    z = (z+nz) % nz;
                    indin = z*nxyin + iin*(nx/2+1) + jin;
                    val1re = band1re[indin];
                    if (order1>0)
                      val1im = band1im[indin];
                    else
                      val1im = 0; // just to suppress a compiling warning
                    root = sqrt(mag2(otf1) + mag2(otf12));/* this scaling is correct if we were to do the linear */
                    fact = otf12/root;                      /*  regression here in fourier space. Guess it is */
                                           /*  right also after transforming to real space  */

                    if (conj) {
                      val1re = conjf(val1re);                    /* do conjugate before multiplying fact */
                      if (order1>0)
                        val1im = conjf(val1im);
                    }
                    val1re *= fact;
                    if (order1>0)
                      val1im *= fact;

                    z = (z0+nz) % nz;
                    indout = z*nxy + i*nx + j;
                    if (order1>0)
                      overlap0[indout] = val1re + I*val1im;
                    else
                      overlap0[indout] = val1re;
                  }
                }
              }
            }
          }
          if (rdist21<=rdistcutoff) { /* within theoretical OTF support of the center band, if x1,y1 are in the side band frame */

            for (z0=-zdistcutoff;z0<=zdistcutoff ;z0++) { /* Doing the loop in a backwards z-fastest order to avoid */
              /* repeating the dist calculations, and the empty xy corners. Is this stupid? */
              /* Is the caching advantage of the natural x-fastest order more important? */
              if (!(z0==0 && pParams->bNoKz0)) { /* use kz==0 plane or not? */
                otf2 = otfinterpolate(OTF[order2], x1, y1, z0, kzscale, pParams);
                if (cabsf(otf2)*order0_2_factor>otfcutoff) {  /* within the area of usably large OTF values of the side band that is to be scaled */

                  otf21 = otfinterpolate(OTF[order1], x21, y21, z0, kzscale, pParams);
                  if (cabsf(otf21)>otfcutoff) { /* within usable OTF support of the center band, if x1,y1 are in the side band frame */
                    if (conj) z = -z0;
                    else z = z0;
                    z = (z+nz) % nz;
                    indin = z*nxyin + iin*(nx/2+1) + jin;
                    val2re = band2re[indin];
                    val2im = band2im[indin];
                    root = sqrt(mag2(otf2) + mag2(otf21)); /* this scaling is correct if we were to do linear */
                    fact = otf21/root;                      /*  regression here in fourier space. Guess it is */
                                          /*  right also after transforming to real space  */

                    if (conj) {
                      val2re = conjf(val2re);
                      val2im = conjf(val2im);
                    }                      /* do conjugate before multiplying fact */
                    val2re *= fact;
                    val2im *= fact;

                    z = (z0+nz) % nz;
                    indout = z*nxy + i*nx + j;
                    overlap1[indout] = val2re + I*val2im;  /* bandplus = bandre + i*bandim */
                  }
                }
              }
            }
          }
        }
      }
    }  /* done making the arrays */


  fftwf_execute(fftplan);     /*** Must use full complex fft ***/
  fftwf_execute_dft(fftplan, overlap1, overlap1);
  /* arrays are now in real space, but are still complex since we cut out a non-hermitian part in fourier space */
  fftwf_destroy_plan(fftplan);
}



/************************* fitk0andmodamps *******************************/
/*     Refines the k0-vector (wave vector of the illumination structure) */
/*     and finds the modulation amplitude and phase in real space by     */
/*     iterating a complex linear regression of bandplus vs. band0       */
/*     (after having transformed to real space only the overlapping      */
/*     parts of the bands). First adjusts the angle of k0, assuming      */
/*     its magnitude is well known, and then refines  the magnitude.     */
/*************************************************************************/

void fitk0andmodamps(fftwf_complex *bands[], fftwf_complex *overlap0, fftwf_complex *overlap1, int nx, int ny, int nz, int norders, vector *k0, float dxy, float dz, fftwf_complex *otf[], short wave, fftwf_complex amps[], ReconParams * pParams)
{
  float k0mag, k0angle, mag, angle, amp1=0., amp2, amp3, x1=0., x2, x3, a;
  float deltaangle = 0.001, deltamag;  //deltamag in 1/micron
  fftwf_complex modamp;
  int redoarrays, fitorder1, fitorder2, order;

  fitorder1 = 0;
  if (nz>1) {
    if (!pParams->bBessel)
      fitorder2 = 2;
    else
      fitorder2 = 1;
  }
  else
    fitorder2 = 1;

  k0mag = sqrt(k0->x * k0->x + k0->y * k0->y);  // in 1/micron
  k0angle = atan2(k0->y, k0->x);

  deltamag = 0.1 / ( (nx>ny ? nx: ny) * dxy );

  redoarrays = (pParams->recalcarrays>=1);    /* recalculate the overlap arrays at least this first time */
  x2 = k0angle;    amp2 = getmodamp(k0angle, k0mag, bands, overlap0,  overlap1, nx, ny, nz, fitorder1, fitorder2, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 0);

  redoarrays = (pParams->recalcarrays>=3);    /* recalculate the overlap arrays every time only if recalcarrays >= 3*/
  angle = k0angle + deltaangle;
  x3 = angle;      amp3 = getmodamp(angle, k0mag, bands, overlap0,  overlap1, nx, ny, nz, fitorder1, fitorder2, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
  if (amp3>amp2) {
    while(amp3>amp2) {
      amp1 = amp2;  x1 = x2;
      amp2 = amp3;  x2 = x3;
      angle += deltaangle;
      x3 = angle;      amp3 = getmodamp(angle, k0mag, bands, overlap0, overlap1, nx, ny, nz, fitorder1, fitorder2, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
    }
  }
  else {
    angle = k0angle;
    a=amp3; amp3=amp2; amp2=a;
    a=x3; x3=x2; x2=a;
    while(amp3>amp2) {
      amp1 = amp2;  x1 = x2;
      amp2 = amp3;  x2 = x3;
      angle -= deltaangle;
      x3 = angle;      amp3 = getmodamp(angle, k0mag, bands, overlap0, overlap1, nx, ny, nz, fitorder1, fitorder2, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
    }
  }  /* the maximum of modamp(x) is now between x1 and x3 */
  angle = fitxyparabola(x1, amp1, x2, amp2, x3, amp3);   /* this should be a good angle.  */

  /***** now search for optimum magnitude, at this angle  *****/

  x2 = k0mag;    amp2 = getmodamp(angle, k0mag, bands, overlap0, overlap1, nx, ny, nz, fitorder1, fitorder2, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 0);

  mag = k0mag + deltamag;
  x3 = mag;      amp3 = getmodamp(angle, mag, bands, overlap0, overlap1, nx, ny, nz, fitorder1, fitorder2, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
  if (amp3>amp2) {
    while(amp3>amp2) {
      amp1 = amp2;  x1 = x2;
      amp2 = amp3;  x2 = x3;
      mag += deltamag;
      x3 = mag;      amp3 = getmodamp(angle, mag, bands, overlap0, overlap1, nx, ny, nz, fitorder1, fitorder2, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
    }
  }
  else {
    mag = k0mag;
    a=amp3; amp3=amp2; amp2=a;
    a=x3; x3=x2; x2=a;
    while(amp3>amp2) {
      amp1 = amp2;  x1 = x2;
      amp2 = amp3;  x2 = x3;
      mag -= deltamag;
      x3 = mag;      amp3 = getmodamp(angle, mag, bands, overlap0, overlap1, nx, ny, nz, fitorder1, fitorder2, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
    }
  }  /* the maximum of modamp(x) is now between x1 and x3 */

  mag = fitxyparabola(x1, amp1, x2, amp2, x3, amp3);  /* this should be a good magnitude.  */
  /* if we were perfectionist we'd iterate for angle again now */

  printf("Optimum modulation amplitude:\n");
  redoarrays = (pParams->recalcarrays>=2);    /* recalculate the overlap arrays for optimum modamp fit */
  amp3 = getmodamp(angle, mag, bands, overlap0,  overlap1, nx, ny, nz, fitorder1, fitorder2, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 1);
  /* one last time, to find the modamp at the optimum k0*/

  printf("Optimum k0 angle=%f, length=%f (1/microns), spacing=%f microns\n", angle, mag, 1.0/mag);

  k0->x = mag*cos(angle); k0->y = mag*sin(angle);
  amps[fitorder2] = modamp;

  /* finally find the modamp for the other orders */
  redoarrays=1;
  if (nz==1)
    for (order=2; order<norders; order++) {
      /* assuming that "angle" and "mag" remain the same for every adjacent pair of bands within one direction */
      getmodamp(angle, mag, bands, overlap0, overlap1, nx, ny, nz, order-1, order, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 1);
      amps[order] = modamp;
    }
  else /* 3D */
    for (order=1; order<norders; order++)
      if (order!=fitorder2) {
        getmodamp(angle, mag, bands, overlap0, overlap1, nx, ny, nz, 0, order, dxy, dz, otf, wave, &modamp, redoarrays, pParams, 1);
        amps[order] = modamp;
      }

  return;
}



/***************************** getmodamp ************************************/
/*     Allows findrealspacemodamp to be called with angle and magnitude     */
/*     of the k0 vector instead of its x and y components.  Returns as      */
/*     the function value the squared magnitude of the fitted modulation    */
/*     amplitude, but also yields its full complex value as the second last */
/*     parameter.                                                           */
/****************************************************************************/
float getmodamp(float kangle, float klength, fftwf_complex *bands[], fftwf_complex *overlap0, fftwf_complex *overlap1, int nx, int ny,int nz, int order1, int order2, float dxy, float dz, fftwf_complex *otf[], short wave, fftwf_complex *modamp, int redoarrays, ReconParams *pParams, int bShowDetail)
{
  vector k1;
  float amp2, corr_coef;
  fftwf_complex amp_inv, amp_combo;

  k1.x = klength*cos(kangle); k1.y = klength*sin(kangle);
  corr_coef = findrealspacemodamp(bands, overlap0, overlap1, nx, ny, nz, order1, order2, k1, dxy, dz, otf,
                                  wave, modamp, &amp_inv, &amp_combo, redoarrays, pParams);
  amp2 = mag2(*modamp);

  printf(" In getmodamp: angle=%f, mag=%f, amp=%f, phase=%f\n",kangle, klength, sqrt(amp2), cargf(*modamp));
  if (bShowDetail) {
    printf(" Reverse modamp is: amp=%f, phase=%f\n", 1/cabsf(amp_inv), -cargf(amp_inv));
    printf(" Combined modamp is: amp=%f, phase=%f\n", cabsf(amp_combo), cargf(amp_combo));
    printf(" Correlation coefficient is: %f\n", corr_coef);
  }

  return(amp2);
}



/***************************** makemodeldata **********************************/
/*     Makes model data by constructing a perfect point source corresponding  */
/*     to the OTF.                                                            */
/* Note: this only works with even number dimensions, I think...              */
/******************************************************************************/

void makemodeldata(int nx, int ny, int nz, fftwf_complex *bands[], int norders, vector k0, float dkx, float dky, float dz, fftwf_complex *OTF[], short wave, ReconParams *pParams)
{
  float kzscale, dist1, dkz, rdistcutoff;
  fftwf_complex otf1, *bandptr=0, * bandptr2=0;
  int i, j, k, x1, y1, z1, zdistcutoff[MAXPHASES/2+1];
  float phi[MAXPHASES/2+1]={0.0,0.,0.,1.2,1.8}, sign, modamp[MAXPHASES/2+1]={1.0,1.0,1.0,0.4,0.2};
  int order, indout, nxy;
  float lambdaem, lambdaexc, beta, k0mag;

  if (dz>0)
    dkz = 1/(nz*dz);
  else
    dkz = pParams->dkzotf;
  kzscale = dkz / pParams->dkzotf;
  k0mag = sqrt(k0.x*k0.x + k0.y*k0.y);   /* k0 magnitude (for highest order) in inverse microns */
  lambdaem = (wave/pParams->nimm)/1000.0;  /* emission wavelength in the sample, in microns */
  lambdaexc = 0.88* lambdaem;;  /* 0.88 approximates a typical lambdaexc/lambdaem  */
  beta = asin(k0mag/(2/lambdaexc));   /* angle of side illumination beams */
  rdistcutoff = pParams->na*2/(wave/1000.0);    /* OTF support limit in data pixels */

  if (rdistcutoff> nx/2*dkx) rdistcutoff = nx/2*dkx;

  /* 080201: zdistcutoff[0] depends on options -- single or double lenses */
  if (!pParams->bTwolens) {
    zdistcutoff[0] = (int)ceil(((1-cos(asin(pParams->na/pParams->nimm)))/((wave/pParams->nimm)/1000.0)) / dkz);    /* OTF support axial limit in data pixels */
    zdistcutoff[norders-1] = 1.3*zdistcutoff[0];    /* approx max axial support limit of the OTF of the high frequency side band */
    if (norders>=3)
      for (order=1;order<norders-1;order++)
        zdistcutoff[order] = 2*zdistcutoff[0];       /* axial support limit of the OTF of the medium frequency side band(s?) */
  }
  else {
    zdistcutoff[0] = (int) ceil((2/lambdaem + 2/lambdaexc) / dkz);
    zdistcutoff[norders-1] = (int) ceil((2/lambdaem + 2*cos(beta)/lambdaexc) / dkz);    /* approx max axial support limit of the OTF of the high frequency side band */
    if (norders>=3)
      for (order=1;order<norders-1;order++) {
        float a;
        a = ((float)order)/(norders-1);
        zdistcutoff[order] = (1-a)*zdistcutoff[0] + a*zdistcutoff[norders-1];       /* axial support limit of the OTF of the medium frequency side band(s?) */
      }
  }


  for (order=0;order<norders;order++)
    if (zdistcutoff[order]>=nz/2) zdistcutoff[order]=((nz/2-1) > 0 ? (nz/2-1) : 0);

  nxy = (nx/2+1)*ny;

  for (i=0; i<2*norders-1; i++)
    memset(bands[i], 0, nxy*nz*sizeof(fftwf_complex));
  for (order=0; order<norders; order++)
    for (k=-zdistcutoff[order]; k<=zdistcutoff[order]; k++)
      for (i=-ny/2;i<ny/2;i++)
        for (j=0;j<nx/2+1;j++) {
          x1 = j; y1 = i; z1 = k;
          if (y1<0) y1+=ny;
          if (z1<0) z1+=nz;

          dist1 = sqrt(i*dky*i*dky+j*dkx*j*dkx);
          indout = z1*nxy+y1*(nx/2+1)+x1;

          if (order==0) bandptr = bands[0];
          else {
            bandptr = bands[2*order-1];
            bandptr2 = bands[2*order];
          }

          if (dist1>rdistcutoff) {
            bandptr[indout] = 0;
            if (order >0)
              bandptr2[indout] = 0;
          }
          else {
            sign = 1.0;
            if ((x1+y1+z1)%2 != 0) sign = -1.0;
            otf1 = otfinterpolate(OTF[order], j*dkx, i*dky, k, kzscale, pParams);
            if (order ==0)
              bandptr[indout] = otf1 * sign;
            else {
              bandptr[indout]  = otf1 * cos(phi[order]) *sign *modamp[order];
              bandptr2[indout] = otf1 * sin(phi[order]) *sign *modamp[order];
            }
          }
        }

  return;
}



/********************* findrealspacemodamp  *******************************/
/*     Finds the modulation amplitude and phase in real space by          */
/*     complex linear regression of bandplus vs. band0, after having      */
/*     transformed to real space only the overlapping parts of the bands. */
/*     070699:  Uses true floating point  k0 vector.                      */
/**************************************************************************/

float findrealspacemodamp(fftwf_complex *bands[], fftwf_complex *overlap0, fftwf_complex *overlap1, int nx, int ny, int nz, int order1, int order2, vector k0, float dxy, float dz, fftwf_complex *OTF[], short wave, fftwf_complex *modamp1, fftwf_complex *modamp2, fftwf_complex *modamp3, int redoarrays, ReconParams *pParams)
{
  fftwf_complex amp, Xval, Yval, *expiphi;
  int nxy,z, i,j,xcent,ycent,indxy,ind;
  float angle, kx, ky;
  double sumXmag, sumYmag;
  complex double sumXstarY, sumYstarX;  /* double precision complex */

  complex double sumX, sumY, sumXY;
  fftwf_complex Xmean, Ymean;
  float corr_coef, sumX2, sumY2, tan2beta, beta, modamp_amp, modamp_arg;

  nxy=nx*ny;
  if (redoarrays)
    /* make arrays that contain only the overlapping parts of fourier space.
       Otf-equalize there, set to zero elsewhere  */
    makeoverlaps(bands, overlap0, overlap1, nx, ny, nz, order1, order2, k0.x, k0.y, dxy, dz ,OTF, wave, pParams);

  /* arrays are now in real space, but are still complex since we cut out a non-hermitian part in fourier space */

  xcent = nx/2;  ycent = ny/2;
  sumXstarY = 0.0; sumXmag = 0.0;
  sumYstarX = 0.0; sumYmag=0;

  /* Multiplying k0 by (order2 - order1) is necessary for the 3D case where order2=2 and order1=0; */
  /* in 2D cases, order2 is almost always order1+1, so multiplying by (order2-order1) won't matter. */
  kx = k0.x * (order2 - order1);
  ky = k0.y * (order2 - order1);
  expiphi = (fftwf_complex *) malloc(nx*ny*sizeof(fftwf_complex));
  if (!expiphi) {
    printf("Cannot allocate memory in findrealspacemodamp(). Program aborted\n");
    exit(-1);
  }

#pragma omp parallel for private(i, j, angle) shared(expiphi)
  for (i=0; i<ny; i++)
    for (j=0; j<nx; j++) {
      /* angle = 2*M_PI * ((j-xcent)*kx/nx + (i-ycent)*ky/ny); */
      angle = 2*M_PI * ((j-xcent)*kx + (i-ycent)*ky) * dxy;
      expiphi[i*nx+j] = cos(angle) + I * sin(angle);
    }

// Try reduction in OpenMP for this loop:
#pragma omp parallel for default(shared) reduction(+:sumXstarY, sumXmag, sumYstarX, sumYmag) private(z, i, j, indxy, ind, Xval, Yval)
  for (z=0; z<nz; z++) {
    int indz = z*nxy;

    for (i=0;i<ny;i++)
      for (j=0;j<nx;j++) {
        indxy = i*nx+j;
        ind = indz + indxy;
        Xval = overlap0[ind];
        Yval = overlap1[ind] * expiphi[indxy];

        sumXstarY += conjf(Xval) * Yval;

        sumXmag += mag2(Xval);

        sumYstarX += Xval * conjf(Yval);
        sumYmag += mag2(Yval);
      }
  }

  amp = sumXstarY / sumXmag;

  *modamp1 = amp;

  amp = sumYstarX / sumYmag;

  *modamp2 = amp;

  /* Mats's new measure of modamp */
  tan2beta = 2*cabs(sumXstarY)/(sumXmag - sumYmag);
  beta = 0.5*atan(tan2beta);
  if (beta < 0) beta += M_PI/2;
  modamp_amp = tan(beta);
  modamp_arg = carg(sumXstarY);
  *modamp3 = modamp_amp * (cos(modamp_arg) + I * sin(modamp_arg));

  /* The following is the calculation of correlation coefficients between X and Y */
  /* this loop is to find out the mean of X and Y */
  sumX = 0;
  sumY = 0;
  sumXY = 0;
  sumX2 = 0;
  sumY2 = 0;

#pragma omp parallel for default(shared) reduction(+:sumX, sumY) private(z, i, j, indxy, ind, Xval, Yval)
  for (z=0; z<nz; z++) {
    int indz = z*nxy;

    for (i=0;i<ny;i++)
      for (j=0;j<nx;j++)  {
        indxy = i*nx+j;
        ind = indz + indxy;
        Xval = overlap0[ind];
        Yval = overlap1[ind] * expiphi[indxy];

        sumX += Xval;
        sumY += Yval;
      }
  }
  Xmean = sumX / (nz*nxy);
  Ymean = sumY / (nz*nxy);

  /* the loop is to get the terms needed for calculating the correlation coefficient */
  /* Lin 08/01/2009: looks like Xmean and Ymean are calculated but not used */
#pragma omp parallel for default(shared) reduction(+:sumXY, sumX2, sumY2) private(z, i, j, indxy, ind, Xval, Yval)
  for (z=0; z<nz; z++) {
    int indz = z*nxy;
    fftwf_complex X1, Y1;

    for (i=0;i<ny;i++)
      for (j=0;j<nx;j++)  {
        indxy = i*nx+j;
        ind = indz + indxy;
        Xval = overlap0[ind];
        Yval = overlap1[ind] * expiphi[indxy];

        /*      X1.re = Xval.re - Xmean.re; */
        /*      X1.im = Xval.im - Xmean.im; */
        /*      Y1.re = Yval.re - Ymean.re; */
        /*      Y1.im = Yval.im - Ymean.im; */
        X1=Xval;
        Y1=Yval;
        sumXY += conjf(X1) * Y1;
        sumX2 += mag2(X1);
        sumY2 += mag2(Y1);
      }
  }

  corr_coef = mag2(sumXY) / (sumX2 * sumY2);

  free(expiphi);
  return(sqrt(corr_coef));

}

float order0damping(float radius, float zindex, float rlimit, int zlimit)
{
  float rfraction, zfraction;

  rfraction = radius/rlimit;
  zfraction = fabs(zindex/zlimit);

  return rfraction*rfraction + zfraction*zfraction*zfraction;
}


/*****************************  filterbands  ******************************/
/*     Filters the bands as if they were to be assembled into a (larger)  */
/*     array in fourier space, but does not actually assemble them.       */
/*     (Assembly will take place later in real space.)  Simultaneously    */
/*     applies a wiener-type filter so that the array is ready to be      */
/*     transformed back to real space.                                    */
/**************************************************************************/
void filterbands(int dir, fftwf_complex *bands[], vector k0[], int ndirs, int norders, fftwf_complex *OTF[], float dxy, float dz, fftwf_complex *amp[], float *noiseVarFactors, int nx, int ny, int nz, short wave, ReconParams* pParams)
{
  fftwf_complex *conjamp, *bandptr=0, *bandptr2=0;
  float rdistcutoff, *ampmag2;
  int order, x1, y1, *zdistcutoff, xcent, ycent, nxy, nxyin;
  float apocutoff, zapocutoff, kx, ky, dkx, dky, dkz, kzscale, k0mag;
  float lambdaem, lambdaexc, alpha, beta, betamin, wiener;

  float min_dkr, suppRadius;
  fftwf_complex *tempbandplus = calloc(nx*ny*nz, sizeof(fftwf_complex));

  wiener = pParams->wiener*pParams->wiener;
  xcent = nx/2;  ycent = ny/2;
  nxy = nx*ny;
  nxyin = (nx/2+1)*ny;
  /* dkr = (1/(ny*dy));   /\* inverse microns per pixel in data *\/ */
  if (dz>0)
    dkz = (1/(nz*dz));   /* inverse microns per pixel in data */
  else
    dkz = pParams->dkzotf;
  /* krscale = dkr / pParams->dkrotf;   /\* ratio of radial direction pixel scales of data and otf *\/ */
  kzscale = dkz / pParams->dkzotf;   /* ratio of axial direction pixel scales of data and otf */
  /* k0pix =  sqrt(k0[0].x*k0[0].x + k0[0].y*k0[0].y);   /\* k0 magnitude (for highest order) in pixels *\/ */
  k0mag = sqrt(k0[0].x*k0[0].x + k0[0].y*k0[0].y);   /* k0 magnitude (for highest order) in inverse microns */
  lambdaem = (wave/pParams->nimm)/1000.0;  /* emission wavelength in the sample, in microns */
  lambdaexc = 0.88* lambdaem;;  /* 0.88 approximates a typical lambdaexc/lambdaem  */
  alpha = asin(pParams->na/pParams->nimm);  /* aperture angle of objectives */
  beta = asin(k0mag/(2/lambdaexc));   /* angle of center of side illumination beams */
  betamin = asin((k0mag/(2/lambdaexc)) -sin(alpha)*SPOTRATIO);   /* angle of inner edge of side illumination beams */
  rdistcutoff = (pParams->na*2/(wave/1000.0));    /* OTF support radial limit in 1/micron */
  if (rdistcutoff > 1/(2*dxy)) rdistcutoff = 1/(2*dxy);

  /* 080201: zdistcutoff[0] depends on options -- single or double lenses */
  zdistcutoff = (int *) malloc(norders * sizeof(int));
  if (!pParams->bTwolens && !pParams->bBessel) {
    zdistcutoff[0] = (int) ceil(((1-cos(alpha))/lambdaem) / dkz);    /* OTF support axial limit in data pixels */
    zdistcutoff[norders-1] = 1.3*zdistcutoff[0];    /* approx max axial support limit of the OTF of the high frequency side band */
    if (norders>=3)
      for (order=1;order<norders-1;order++)
        zdistcutoff[order] = (1+lambdaem/lambdaexc)*zdistcutoff[0];       /* axial support limit of the OTF of the medium frequency side band(s?) */
  }
  else if (pParams->bBessel) {
    float kzExMax, halfangle;
    kzExMax = 2 *BesselNA / BesselLambdaEx;
    zdistcutoff[0] = (int) ceil((kzExMax + (1-cos(alpha))/lambdaem) / dkz);    /* OTF support axial limit in data pixels */
    for (order=1;order<norders-1;order++) { 
      halfangle = acos(k0mag*order / kzExMax);
      zdistcutoff[order] = (int) ceil((kzExMax * sin(halfangle) + (1-cos(alpha))/lambdaem) / dkz);
    }
  }
  else {  /* two lenses */
    zdistcutoff[0] = (int) ceil(1.02*(2/lambdaem + 2/lambdaexc) / dkz);  /* 1.02 is just a safety margin */
    zdistcutoff[norders-1] = (int) ceil(1.02*(2/lambdaem + 2*cos(beta)/lambdaexc) / dkz);    /* approx max axial support limit of the OTF of the high frequency side band */
    if (norders==3) {
      zdistcutoff[1] =  (int) ceil(1.02*(2/lambdaem + (1+cos(betamin))/lambdaexc) / dkz); /* axial support limit of the OTF of the medium frequency side band */
    }
    else if (norders>3)
      for (order=1;order<norders-1;order++) {
        float a;
        a = ((float)order)/(norders-1);
        zdistcutoff[order] = 1.1*((1-a)*zdistcutoff[0] + a*zdistcutoff[norders-1]);       /* axial support limit of the OTF of the medium frequency side bands */ /* 1.1 is a blur margin */
      }
  }

  for (order=0;order<norders;order++) {
    if (zdistcutoff[order]>=nz/2) zdistcutoff[order]=((nz/2-1) > 0 ? (nz/2-1) : 0);
    /* printf("order=%d, rdistcutoff=%f, zdistcutoff=%d\n", order, rdistcutoff, zdistcutoff[order]); */
  }

  apocutoff = rdistcutoff + k0mag * (norders-1);

  if (pParams->bTwolens || pParams->bBessel)
    zapocutoff = zdistcutoff[0];
  else
    zapocutoff = zdistcutoff[1];

  ampmag2 = (float *) malloc(norders * sizeof(float));
  conjamp = (fftwf_complex *) malloc(norders * sizeof(fftwf_complex));
  for (order=1;order<norders;order++) {
    ampmag2[order] = mag2(amp[dir][order]);
    conjamp[order] = conjf(amp[dir][order]);
  }

  dkx = 1./(nx*dxy);
  dky = 1./(ny*dxy);
  min_dkr = dkx < dky ? dkx : dky;
  suppRadius = pParams->suppression_radius*min_dkr;

  for (order=0;order<norders;order++) {
    if (order==0) bandptr = bands[0];
    else {
      bandptr = bands[2*order-1];     /* bands contains only the data of one direction -- dir*/
      bandptr2 = bands[2*order];
    }
    kx = order * k0[dir].x;
    ky = order * k0[dir].y;

#pragma omp parallel for private(y1, x1)
    for (y1=-(ycent-1);y1<=ycent;y1++) { /* (integer) coords of an imagined full four-quadrant array of the band to be filtered, */
      float weight, sumweight, amp2mag2, rdist1, rdist2;
      fftwf_complex bandreval, bandimval, bandplusval, otf1, otf2, scale;
      float xabs, yabs, x2, y2, kx2, ky2, rdistabs, dampfact, apofact=1.;
      int order2, dir2, iin, jin, conj, xyind, ind, z0, iz, z;
      int iout, jout, zout; // for tempbandplus indexing

      for (x1=-(xcent-1);x1<=xcent;x1++) { /* with origin of fourier space at (0,0) of each band*/
        float x1f, y1f;   // (x1, y1) counterpart in 1/micron
        /*x1, y1 are coords within each band to be scaled */
        if (x1>=0) {
          iin = y1;   /* integer coords of actual arrays to be filtered */
          jin = x1;
          conj = 0;
        } else {
          iin = -y1;
          jin = -x1;
          conj = 1;
        }
        if (order==0 && conj) continue;   /* ?? For center band only the actual pixels need to be filtered?  */
        if (iin<0) iin+=ny;
        xyind = iin*(nx/2+1)+jin;

        x1f = x1 * dkx;
        y1f = y1 * dky;
        rdist1 = sqrt(x1f*x1f+y1f*y1f);  /* dist from center of band to be filtered */
        if (rdist1<=rdistcutoff) {   /* is x1,y1 within the theoretical lateral OTF support of the data that is to be scaled? */
          xabs = x1f + kx;   /* coords (in 1/micron) relative to absolute Fourier space, with */
          yabs = y1f + ky;   /* the absolute origin=(0,0) after the band is shifted by k0 */
          rdistabs = sqrt(xabs*xabs + yabs*yabs);  // used later for apodization calculation
          for (z0=-zdistcutoff[order];z0<=zdistcutoff[order];z0++) {
            otf1 = otfinterpolate(OTF[order], x1f, y1f, z0, kzscale, pParams);
            weight = mag2(otf1);
            if (order!= 0) weight *= ampmag2[order];
            dampfact = 1. / noiseVarFactors[dir*norders+order];
            if (pParams->bSuppress_singularities && order != 0 && rdist1 <= suppRadius)
              dampfact *= suppress(rdist1, suppRadius);

            else if (!pParams->bDampenOrder0 && !pParams->bNoOrder0 && pParams->bSuppress_singularities
                     && order ==0 && rdist1 <= suppRadius)
              dampfact *= suppress(rdist1, suppRadius);
                    
            else if (pParams->bDampenOrder0 && order ==0)
              dampfact *= order0damping(rdist1, z0, rdistcutoff, zdistcutoff[0]);

            else if (pParams->bNoOrder0 && order ==0)
              dampfact = 0.;

            // if no kz=0 plane is used:
            if (order==0 && z0==0 && pParams->bNoKz0) dampfact = 0;

            weight *= dampfact;
            sumweight=weight;

            for (dir2=0; dir2<ndirs; dir2++) {
              for (order2=-(norders-1);order2<norders;order2++) {
                if (dir2==dir && order2==order) continue;
                if (!pParams->bFilteroverlaps && !(order2==0 && order==0)) continue; /* filteroverlaps is always true except when (during debug) generating an unfiltered exploded view */
                amp2mag2 = mag2(amp[dir2][abs(order2)]);
                kx2 = order2 * k0[dir2].x;
                ky2 = order2 * k0[dir2].y;
                x2 = xabs-kx2; /* coords rel to shifted center of band 2 */
                y2 = yabs-ky2;
                rdist2 = sqrt(x2*x2 + y2*y2);       /* dist from center of band 2 */
                if (rdist2<rdistcutoff) {
                  otf2 = otfinterpolate(OTF[abs(order2)], x2, y2, z0, kzscale, pParams);
                  weight = mag2(otf2) / noiseVarFactors[dir2*norders+abs(order2)];
                  if (order2 != 0) weight *= amp2mag2;

                  if (pParams->bSuppress_singularities && order2 != 0 && rdist2 <= suppRadius)
                    weight *= suppress(rdist2, suppRadius);

                  else if (!pParams->bDampenOrder0 && !pParams->bNoOrder0 && pParams->bSuppress_singularities
                           && order2 ==0 && rdist2 <= suppRadius)
                    weight *= suppress(rdist2, suppRadius);
                    
                  else if (pParams->bDampenOrder0 && order2 ==0)
                    weight *= order0damping(rdist2, z0, rdistcutoff, zdistcutoff[0]);
                  else if (pParams->bNoOrder0 && order2 ==0)
                    weight = 0.;

                  if (pParams->bNoKz0 && order2==0 && z0==0) weight = 0;

                  sumweight += weight;
                  
                }
              }
            }
            sumweight += wiener;
            scale = dampfact * conjf(otf1) / sumweight;

            if (pParams->apodizeoutput) {
              float rho, zdistabs;
              zdistabs = abs(z0);

              if (zapocutoff > 0) {  /* 3D case */
                if (!pParams->bBessel)
                  rho = sqrt((rdistabs/apocutoff)*(rdistabs/apocutoff)+(zdistabs/zapocutoff)*(zdistabs/zapocutoff));
                else {
                  float rhox, rhoy, rhoz;
                  rhox = xabs/apocutoff * 8./5.;
                  rhoy = yabs/apocutoff;
                  rhoz = zdistabs/zapocutoff;
                  rho = sqrt(rhox*rhox + rhoy*rhoy + rhoz*rhoz);
                }
              }
              else         /* 2D case */
                if (pParams->bPreciseApoBoundary) {
                  float pointArg = atan2(yabs, xabs);
                  rho = rdistabs/getPreciseResLimit(pointArg, boundaryAngles, k0,
                                                    ndirs, norders, rdistcutoff);
                }
                else
                  rho = rdistabs/apocutoff;

              if (rho > 1) rho = 1;

              if (pParams->apodizeoutput == 1)    /* cosine-apodize */
                apofact = cos((M_PI*0.5)* rho);
              else if (pParams->apodizeoutput == 2) /* apodizing by a power function */
                apofact = pow(1.0 - rho, pParams->apoGamma);
              scale *= apofact;
            }

            /* What we want is to use mag2 for the weights, as you have done, and then
               set  scale = conjugate(otf1)/sumweight */

            /* separate (for this pixel) the even and odd "bands" into the true plus and minus bands */
            /* apply the scale to the plus band only (this will apply it to the minus band as well by symmetry?) */
            /* reassemble into the even and odd "bands" */
            if (conj) z = -z0;
            else z = z0;
            iz = (z+nz)%nz;   /* coords of the fourier space arrays have the origin of fourier space at (0,0,0) */
            ind = iz*nxyin+xyind;   /* index of the corresponding point in the input array */
            if (order == 0) {
              if (!conj)
                bandptr[ind] *= scale;
              else
                printf("error: order=0 and conj\n");
            }
            else {
                scale *= conjamp[order]; /* not invamp: the 1/|amp| factor is */
                /*  taken care of by including ampmag2 in the weights */
                bandreval = bandptr[ind];
                bandimval = bandptr2[ind];
                if (conj) {
                  bandreval = conjf(bandreval);
                  bandimval = conjf(bandimval);
                }
                bandplusval  = bandreval + I * bandimval;  /* bandplus = bandre + i bandim */
//                bandminusval = bandreval - I * bandimval;  /* bandminus = bandre - i bandim */

                bandplusval *= scale;   /* scale only the bandplus part - bandminus will take care */
                /* of itself because of symmetry (?) */

                iout = (y1+ny)%ny;
                jout = (x1+nx)%nx;
                zout = (z0+nz)%nz;
                tempbandplus[zout*nxy+iout*nx+jout] = bandplusval;

                /* bandreval = 0.5*(bandplusval + bandminusval); */
                /* bandimval = - 0.5*(bandplusval - bandminusval) * I; */
                /* if (conj) { */
                /*   bandreval = conjf(bandreval); */
                /*   bandimval = conjf(bandimval); */
                /* } */
                /* bandptr[ind] = bandreval; */
                /* bandptr2[ind] = bandimval; */
            }
          }   /* for z0=... */
        }
        else { /* if rdist1>rdistcutoff */

          for (z0=-zdistcutoff[order]; z0<=zdistcutoff[order]; z0++) {
            iz = (z0+nz)%nz;   /* coords of the fourier space arrays have the origin of fourier space at (0,0,0) */

            ind = iz*nxyin+xyind;
            bandptr[ind] = 0;
            if (order !=0)
              bandptr2[ind] = 0;
          }
        }
        for (z0=zdistcutoff[order]+1; z0<nz-zdistcutoff[order]; z0++) {   /* zero out all sections between zdistcutoff+1 and
                                                                             nz-zdistcutoff-1 (include) */
          ind = z0*nxyin+xyind;
          bandptr[ind] = 0;
          if (order !=0)
            bandptr2[ind] = 0;
        }
      } // for x1
    }// for y1
    // Convert tempbandplus to bandre, bandim after Muthuvel discussion Feb. 2011
    if(order==0) continue;
#pragma omp parallel for private(y1, x1)
    for (y1=-(ycent-1); y1<=ycent; y1++)
      for (x1=0; x1<=xcent; x1++) {
        int iin, jin, z0, iz, indin, xyind, xyind_full, xyind_full_conj, ind_full, ind_full_conj;
        complex bandplusval, bandminusval;
        iin =(y1+ny) % ny;  // iin, jin index the input (and output) bandptr, and bandptr2
        jin = x1;
        xyind = iin*(nx/2+1)+jin;
        xyind_full = iin*nx + jin; //xyind_full indexes tempbandplus
        xyind_full_conj = (-y1+ny)%ny * nx + (-x1+nx)%nx;
        // xyind_full_conj is the corresponding point of xyind_full in the negative kx half space
        for (z0=-zdistcutoff[order]; z0<=zdistcutoff[order]; z0++) {
          iz = (z0+nz) %nz;
          indin = iz*nxyin+xyind; // indin indexes bandptr, bandptr2
          ind_full = iz*nxy + xyind_full; // ind_full indexes tempbandplus
          ind_full_conj = ((-z0+nz)%nz)*nxy + xyind_full_conj;
          // ind_full_conj is the corresponding point of ind_full in the negative kx half space
          bandplusval = tempbandplus[ind_full];
          bandminusval = conjf(tempbandplus[ind_full_conj]); //bandminusval at ind_full is the conjugate of its value at ind_full_conj
          bandptr[indin] = 0.5*(bandplusval + bandminusval);
          bandptr2[indin] = - 0.5*(bandplusval - bandminusval) * I;
        }
      }
    memset(tempbandplus, 0, nxy*nz*sizeof(fftwf_complex));
  } /* for order */

  free(tempbandplus);
  free(zdistcutoff);
  free(ampmag2);
  free(conjamp);
}


/***************************  Assemblerealspacebands ******************************/
/*  Transforms the bands to real space, multiplies by the right (sub-pixel        */
/*  precision) modulation sine functions, and adds them together.                 */
/*  First data need to be moved to larger arrays, since resolution will increase. */
/**********************************************************************************/
void assemblerealspacebands(int dir, float *outbuffer, fftwf_complex *bigbuffer, fftwf_complex *bands[],
                            int ndirs, int norders, vector k0[], float dxy, int nx, int ny, int nz, 
                            fftwf_plan fftplan, float zoomfact, int z_zoom, float expfact, int bNoOrder0)
{
  int i,j,k,ind,xcent,ycent,indxy;
  float angle,fact,*coslookup,*sinlookup;
  int nxy, order, nxyz;

  coslookup = (float *) malloc((int)(nx*zoomfact*ny*zoomfact*sizeof(float)));
  sinlookup = (float *) malloc((int)(nx*zoomfact*ny*zoomfact*sizeof(float)));
  if (!coslookup || !sinlookup) {
    printf("Cannot allocate memory for sin/cos lookup tables in assemblerealspacebands()\n Program aborted\n");
    exit(-1);
  }

  xcent = (int) ((nx*zoomfact)/2);  ycent = (int)((ny*zoomfact)/2);
  fact = expfact/0.5;  /* expfact is used for "exploded view".  For normal reconstruction expfact = 1.0  */
  nxy = zoomfact*nx*zoomfact*ny;
  nxyz = z_zoom* nz * nxy;

  if (!bNoOrder0) {
    /* move order 0 band to bigbuffer, fill in with zeroes */
    printf("moving centerband\n");
    move(bands, 0, bigbuffer, nx, ny, nz, zoomfact, z_zoom);

    /* transform it */
    printf("re-transforming centerband\n");
    fftwf_execute(fftplan);     /*** full complex fft ***/

    printf("inserting centerband\n");

#pragma omp parallel for private(ind)
    for(ind=0; ind<nxyz; ind++)
      outbuffer[ind] += crealf(bigbuffer[ind]);

    printf("centerband assembly completed\n");
  }
  for (order=1; order < norders; order ++) {
    float k0x, k0y;
    /* move side bands to bigbuffer, fill in with zeroes */
    printf("moving order %d\n",order);
    move(bands, order, bigbuffer, nx, ny, nz, zoomfact, z_zoom);

    /* transform it into real space*/
    fftwf_execute(fftplan);   /*** full complex fft ***/

    /***** For 3D, prepare 2D array of sines and cosines first, then loop over z. ******/
    k0x = k0[dir].x*((float)order);
    k0y = k0[dir].y*((float)order);

#pragma omp parallel for private(i,j,ind,angle) shared(coslookup, sinlookup)
    for (i=0;i<(int)(zoomfact*ny);i++)
      for (j=0;j<(int)(zoomfact*nx);j++) {
        ind = i*zoomfact*nx + j;
        /* angle = fact * M_PI * ((j-xcent)*k0x/(zoomfact*nx) + (i-ycent)*k0y/(zoomfact*ny)); */
        angle = fact * M_PI * ((j-xcent)*k0x + (i-ycent)*k0y) * dxy/zoomfact;
        coslookup[ind] = cos(angle);
        sinlookup[ind] = sin(angle);
      }

#pragma omp parallel for private(k, i, j, indxy, ind) shared(outbuffer, coslookup, sinlookup)
    for (k=0; k<z_zoom*nz; k++) {
      int indz = k*nxy;
      for (i=0; i<zoomfact*ny; i++)
        for (j=0; j<zoomfact*nx; j++) {
          indxy = i*zoomfact*nx + j;
          ind = indz + indxy;
          outbuffer[ind] += crealf(bigbuffer[ind]) * 2*coslookup[indxy] - cimagf(bigbuffer[ind]) * 2*sinlookup[indxy];
        }
    }

    printf("order %d sideband assembly completed\n", order);
  } /* for (order =...) */

  free(coslookup);
  free(sinlookup);

}


/**********************************  Move *******************************************/
/*  Moves the contents of the complex fourier space array "inarray" to the          */
/*  complex fourier space array "outarray."  Inarray has dimensions (nx/2+1)*ny,    */
/*  has the origin of fourier space located at (0,0), and has been                  */
/*  multiplied by a checkerboard of (+1,-1).   Outarray is of dimension             */
/*  (zoomfact*nx)*(zoomfact*ny), has the origin of fourier space at (0,0)           */
/*  and does not have the checkerboard compensation.  Its extra elements            */
/*  get set to zero.                                                                */
/************************************************************************************/
void move(fftwf_complex *inarray[], int order, fftwf_complex *outarray, int nx, int ny, int nz, float zoomfact, int z_zoom)
{
  int     i, j, k, xdim, ydim, zdim, xcent, ycent, zcent, nxy, nxyout, nxyzout;

  xdim=zoomfact*nx; ydim=zoomfact*ny;
  zdim = z_zoom*nz;
  nxy = (nx/2+1)*ny;
  nxyout = xdim*ydim;
  nxyzout = zdim * nxyout;
  xcent = nx/2; ycent = ny/2;
  zcent = nz/2;

#pragma omp parallel for private(k) shared(outarray, nxyzout)
  for (k=0; k<nxyzout; k++)
    outarray[k] = 0;

  /* Mats'version */
  for (k=0; k<nz; k++)   /* kij = (non-centered) coords of an imagined full complex array with dims (nx,ny,nz) and origin of fourier space at (0,0,0) */
#pragma omp parallel for private(i, j) shared(k, outarray, inarray, order, nx, ny, nxy, nxyout, nz, xdim, ydim, zdim, xcent, ycent, zcent)
    for (i=0; i<ny; i++)
      for (j=0; j<nx; j++) {
        int x, y, z, indin, indout, xout, yout, zout, conj;
        fftwf_complex valre, valim, val;
        x = j;   /* xyz = centered coords going between -nx/2+1 to +nx/2 etc */
        if (x>xcent) x-= nx;
        y = i;
        if (y>ycent) y-= ny;
        z = k;
        if (z>zcent) z-= nz;
        xout = x;     /* xout,yot,zout = (non-centered) output coords with zoomed-up dims and origin of fourier space at (0,0,0) */
        if (xout<0) xout += xdim;
        yout = y;
        if (yout<0) yout += ydim;
        zout = z;
        if (zout<0) zout += zdim;
        indout = zout*nxyout + yout*xdim + xout;
        if (x<0) {    /* now xyz get turned into coords of the (half fourier space) input arrays */
          x = -x;
          y = -y;
          z = -z;
          conj = 1;
        }
        else
          conj = 0;

        if (y<0) y += ny;
        if (z<0) z += nz;
        indin = z*nxy + y*(nx/2+1) + x;


        if (order == 0) {
          val = inarray[order][indin];
          if (conj)
            val = conjf(val);
        }
        else {
          valre = inarray[2*order-1][indin];
          valim = inarray[2*order][indin];
          if (conj) {
            valre = conjf(valre);
            valim = conjf(valim);
          }
          val = valre + I * valim;
        }

        outarray[indout] = val;
      }

}
