#include <complex.h>   // Use C99 built-in complex type
#include <fftw3.h>

#define SPOTRATIO 0.1 /*ratio of beam size in pupil to pupil size */

#define MAXPHASES 25

#define K0_WARNING_THRESH 2  /* If k0 initial estimate is off from the guess by over this many pixels, a warning will be displayed */


typedef struct { float x; float y;  } vector;
typedef struct { float x; float y; float z; } vector3d;

typedef struct { double re; double im;  } zomplex; /* double-precision complex */

typedef struct {
  float k0startangle, linespacing;
  int bNegAngleStep;
  float na, nimm;
  int   ndirs, nphases, norders_output;
  float *phaseSteps; /* user-specified non-ideal phase steps, one for each orientation */
  float *phaseList;
  int   bTwolens;    /* whether to process I{^5}S dataset */
  int   bBessel;    /* whether to process Bessel-Sheet SIM dataset */
  int   bFastSIM;   /* fast SIM data is organized differently */

  /* algorithm related parameters */
  float zoomfact;
  int   z_zoom;
  int   nzPadTo;  /* pad zero sections to this number of z sections */
  float explodefact;
  int   bFilteroverlaps;
  int   recalcarrays; /* whether to calculate the overlaping regions between bands just once or always; used in fitk0andmodamps() */
  int   napodize;
  float forceamp[MAXPHASES];
  float *k0angles;
  int   bSearchforvector;
  int   bUseTime0k0;   /* whether to use time 0's k0 fit for the rest of a time series */
  int   apodizeoutput;  /* 0-no apodize; 1-cosine apodize; 2-triangle apodize; used in filterbands() */
  float apoGamma; /* 1.0 means triangular apo; <1.0 means apo curve is above the straight line */
  int   bPreciseApoBoundary; /* whether to precisely define resolution boundary for apodization purpose*/
  int   bSuppress_singularities;  /* whether to dampen the OTF values near band centers; used in filterbands() */
  float   suppression_radius;   /* if suppress_singularities is 1, the range within which suppresion is applied; used in filterbands() */
  int   bDampenOrder0;  /* whether to dampen the OTF values near band centers; used in filterbands() */
  int   bNoOrder0;  /* whether to use order 0 in final assembly; used in filterbands() */
  int   bFitallphases;  /* In NL SIM, whether to use fitted phase for all orders or to infer higher order phase from order 1 phase fil */
  int   do_rescale; /* fading correction method: 0-no correction; 1-with correction */
  int   equalizez;
  int   equalizet;
  int   bNoKz0;   /* if true, no kz=0 plane is used in modamp fit and assemblerealspace() */
  float wiener, wienerInr;
  int   bUseEstimatedWiener;
  int   b2D; /* data is 2D SIM even if nz>1 */

  /* OTF specific parameters */
  int   nxotf, nyotf, nzotf;
  float dkzotf, dkrotf;  /* OTF's pixel size in inverse mirons */
  int   bRadAvgOTF;   /* is radially-averaged OTF used? */

  /* drift correction and phase step correction related flags */
  int   bFixdrift;   /* whether nor not to correct drift between pattern directions */
  float drift_filter_fact; /* fraction of the NA; used as the cutoff frequency in high-pass filtering in drift estimation */
  int   bFix2Ddrift;  /* whether, in addition to fixdrift, to correct drift on a section-by-section basis with each pattern dir*/
  int   bFixPhaseStep;  /* whether to use section-by-section drift estimates for building separation matrix*/
  int   type_of_refImage;  /*Type of reference image to be used in 2D SIM, especially nonlinear SIM --
                           0: use sum of switch-off and exposure images; one per direction
                           1: use phase averaged exposure image; one per direction
                           2: use experimentally measured ref image (switch-off+exposure of the first phase?); one per direction
                           3: use one experimentally measured ref image for the whole dataset (switch-off+exposure of the first phase of the first direction?)
                          */
  int   bFitDrift;    /* Whether to fit the estimated between-exposure 2D drifts with a curve, and use the fitted values in correction */
  int   nOffImages; /* Number of switch-off images in nonlinear SIM; they're taken before every normal exposure images with
                     the grating Pi phase shifted relative to the grating position at exposure; the adidtion of the off images and the
                     normal exposure should be close to a conventional image; i.e., they're "complementary" in some sense.
                  */
  int   combineOffOrders[MAXPHASES];
  int   bCombineOffExp;  /* In NL SIM,  whether to combine switch-off images with the normal exposure images in reconstruction */

  /* Camera related parameters */
  float constbkgd;
  int bBgInExtHdr; /* In Andor EMCCD, background varies with each exposure, esp. in EM mode. Hidden-behind-aluminum-foil pixels can be used to estimate background of each exposure and stored in the extended header. When this option is true, the 3rd float of the extended header stores such estimated background values. */
  int   bUsecorr;    /* whether to use a camera flat-fielding (or correction) file */
  char  corrfiles[200];  /* name of the camera correction file if bUsecorr is 1 */
  float readoutNoiseVar;
  float electrons_per_bit;

  /* Debugging flags */
  int   bMakemodel;  /* whether to fake an ideal point source and obtain a recontruction of it (i.e., an effective PSF in some sense) */
  int   bSaveSeparated; /* whether to save separated bands and then quit before filterbands */
  char  fileSeparated[200];
  int   bSaveAlignedRaw; /* whether to save dirft-corrected raw images (within each direction) and then quit before filterbands */
  char  fileRawAligned[200];
  int   bSaveOverlaps; /* whether to save makeoverlaps() output into a file */
  char  fileOverlaps[200];

  /* performance related parameters */
  int  nthreads;  /* number of threads to use */
} ReconParams;

typedef struct{
  float timestamp, phaseAbsDeg, expDose, xdrift, ydrift;
} myExtHeader;


extern int nthreads;

float mag2(fftwf_complex b);
int image_arithmetic(float * imDest, float *imSrc, int len, float fDest, float fSrc);

void filterbands(int direction, fftwf_complex *bands[], vector *k0, int ndirs, int norders, fftwf_complex *OTF[], float dxy, float dz, fftwf_complex *amp[], float *noiseVars, int nx, int ny, int nz, short wavei, ReconParams * pParams);

void findk0(fftwf_complex * bands[], fftwf_complex *bandplus,fftwf_complex *bandzero,int nx,int ny,int nz,int norders,vector *k0, float dxy, float dz, fftwf_complex *OTF[], short wave, ReconParams * pParams);

fftwf_complex findmodamp(fftwf_complex *bandplus,fftwf_complex *band0,int nx,int ny,int k0x,int k0y, float dy, float *OTF, float dkotf, int notf, short wave);

void  findpeak(float array[], int sizex, int sizey, vector *peak);

float fitparabola( float a1, float a2, float a3 );

float fitxyparabola( float x1, float y1, float x2, float y2, float x3, float y3 );

void fitk0andmodamps(fftwf_complex *bands[], fftwf_complex *overlap0, fftwf_complex *overlap1, int nx, int ny, int nz, int norders, vector *k0, float dxy, float dz, fftwf_complex *otf[], short wave, fftwf_complex amps[], ReconParams *pParams);
                     
float findrealspacemodamp(fftwf_complex *bands[], fftwf_complex *overlap0, fftwf_complex *overlap1, int nx, int ny, int nz, int order1,
                          int order2, vector k0, float dxy, float dz, fftwf_complex *OTF[], short wave,
                          fftwf_complex *modamp1, fftwf_complex *modamp2, fftwf_complex *modamp3, int redoarrays, ReconParams *pParams);


float suppress(float x, float r);

void assemblerealspacebands(int dir, float *outbuffer, fftwf_complex *bigbuffer, fftwf_complex *bands[],
                            int ndirs, int norders, vector k0[], float dxy, int nx, int ny, int nz, fftwf_plan fftplan,
                            float xy_zoomfact, int z_zoomfact, float explodefact, int bNoOrder0);

void move(fftwf_complex **inarray, int order, fftwf_complex *outarray, int nx, int ny, int nz, float zoomfact, int z_zoom);

void rescale(int nx, int ny, int nz, int z, int zoffset, int direction, int wave, int t, int phases, float *floatimage[], int equalizez, int equalizet, double * sum_dir0_phase0);

void makematrix (int nphases, int norders, int dir, float *phases, float *sepMatrix, float * noiseVars);

void apodize(int napodize,int nx,int ny,float *image);

void cosapodize(int nx,int ny,float *image);

void determinedrift_bt_dirs(int ifile, int wave_no, int time_no, vector3d* drift, int nx, int ny, int nz, int ndirs, int nphases, float *background, float dxy, float dz, float rdistcutoff, short wave, ReconParams *pParams);


float getmodamp(float kangle, float klength, fftwf_complex *bands[], fftwf_complex *overlap0, fftwf_complex *overlap1, int nx,int ny, int nz, int order1, int order2, float dxy, float dz, fftwf_complex *otf[], short wave, fftwf_complex *modamp, int redoarrays, ReconParams * pParams, int bShowDetail);

void makeoverlaps(fftwf_complex *bands[], fftwf_complex *overlap0, fftwf_complex *overlap1, int nx, int ny, int nz, int order1, int order2,
                  float k0x, float k0y, float dy, float dz, fftwf_complex *OTF[], short wave, ReconParams *pParams);

void makemodeldata(int nx, int ny, int nz, fftwf_complex *bands[], int norders, vector k0, float dkx, float dky, float dz, fftwf_complex *OTF[], short wave, ReconParams *pParams);


void fixdrift_bt_dirs(fftwf_complex *bands[], int norders, vector3d drift, int nx,int ny, int nz);

void separate(int nx,int ny,int z,int direction,int nphases, int norders, float *floatimage[], float *sepMatrix);

fftwf_complex otfinterpolate(fftwf_complex * otf, float kx, float ky, int kz, float kzscale, ReconParams* pParams);

void load_and_flatfield(int istream_no, int section_no, int wave_no, int time_no, float *bufDestiny, float *buffer, int nx, int ny, float *background, float backgroundExtra, float *slope, float inscale, int bUsecorr);

void getbg_and_slope(char *corrfiles, float *background, float *slope, int nx, int ny);

void SetDefaultParams(ReconParams *pParams);

int commandline(int argc, char *argv[], ReconParams * myParams, char *ifiles, char *ofiles, char * otffiles, int * ifilein, int *ofilein, int * otffilein);

int calcRefImage(float **rawImages, float *refImage, float **offImages, int nOffImages, int nx, int ny, int nphases, int type_of_refImage);

int fitXYdrift(vector3d *drifts, float * timestamps, int nPoints, vector3d *fitted_drift, float *eval_timestamps, int nEvalPoints);

void fixdrift_2D(fftwf_complex **floatimages, vector3d *driftlist, int nphases, int nx, int ny, int nz, int dir, int z);

void determinedrift_2D(fftwf_complex **rawImages, fftwf_complex **offImages, int nOffImages, fftwf_complex * refImage,
                       vector3d *drifts, int nphases, int nx, int ny, int dir,
                       float rdistcutoff, float dkx, float dky, float drift_filter_fact);
void calcPhaseList(float * phaseList, vector3d *driftlist, float *phaseCorrAcq, float k0startangle, float linespacing, float dkr, int nOffImages, int nphases, int nz, int direction, int z);

void combine_exp_off(float **expImages, float **offImages, int nx, int ny, int nphases, int dir,
                     float* noiseVarFactors, float readoutNoiseVar, float electrons_per_bit);
void combine_exp_off_F(fftwf_complex **expImages, fftwf_complex **offImages, int nx, int ny, int norders, int dir, int rdistcutoff, float dkx, float dky, float* noiseVarFactors, int* combineOffOrders);

void matrix_transpose(float * mat, int nRows, int nCols);
void print_matrix( char* desc, int m, int n, float* a, int lda );
float estimate_Wiener(fftwf_complex ** dataF, int nx, int ny, int nz, int nphases, float rdistcutoff);
float estimate_Wiener_3D(fftwf_complex *dataF, int nx, int ny, int nz, float dkx, float dky, float dz, int wavelength, ReconParams* pParams);

double sum_of_pixels(float *im, int nx, int ny, int bExtra2Cols);
void calculateSectorBoundaries(int ndirs, vector* k0, float * boundaryAngles);
float getPreciseResLimit(float arg, float * boundaryAngles, vector *k0, int ndirs, int norders, float R);

void mrc_file_write(const float *buffer, int nx, int ny, int nz, float rlen, float zlen, int mode, int iwave, const char *files);


void determine_otf_dimensions(
#ifdef __SIRECON_USE_TIFF__
                              CImg & otf_tiff,
#else
                              int otfstream_no,
#endif
                              int norders, int nz, ReconParams *pParams, int *sizeOTF);

int loadOTFs(
#ifdef __SIRECON_USE_TIFF__
             CImg & otf_tiff,
#else
             int otfstream_no,
#endif
             fftwf_complex **otf, int norders, int nz, int sizeOTF, ReconParams *pParams);


extern int sgemm_(char *transa, char *transb, int *m, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);

extern int sgetrf_(int *m, int *n, float *a, int *lda, int *ipiv, int *info);

extern int sgetri_(int *n, float *a, int *lda, int *ipiv, float *work, int *lwork, int *info);

extern int sgels_(char *trans, int *m, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, float *work, int *lwork, int *info);
