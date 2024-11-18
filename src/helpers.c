#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <IMInclude.h>

#include "sirecon.h"

int getinteger(char * buffer, int * ret)
{
  int i, len;

  len = strlen(buffer);
  for (i=0; i<len; i++)
	if (!isdigit(buffer[i])) {
	  if (i==0 && buffer[i] !='+' && buffer[i] !='-' )
	    return 0;
	  else if (i!=0)
	    return 0;
	}
  sscanf(buffer, "%d", ret);
  return 1;
}

int getfloat(char *buffer, float *ret)
{
  int i, len;
  int dotcount = 0;

  len = strlen(buffer);
  for (i=0; i<len; i++)
	if (!isdigit(buffer[i])) {
	  if (i==0 && buffer[i] !='+' && buffer[i] !='-' && buffer[i]!='.' )
	    return 0;
	  else if (i!=0 && buffer[i] != '.')
	    return 0;
	  if (buffer[i] == '.')
	    dotcount ++;
	}
  if (dotcount > 1) return 0;
  sscanf(buffer, "%f", ret);
  return 1;

}

int getfloats(char *buffer[], float *ret, int num)
{
  int i;

  for (i=0; i<num; i++) {
	if (!getfloat(buffer[i], ret+i))
	  return 0;
  }
  return 1;
}

int getints(char *buffer[], int *ret, int num)
{
  int i;

  for (i=0; i<num; i++) {
	if (!getinteger(buffer[i], ret+i))
	  return 0;
  }
  return 1;
}

void Usage(char *argv[])
{
  printf("Usage:\n%s [[input file] [output file] [OTF file]] [Options]\n", argv[0]);
  printf("\nOptions:\n");
  printf("\t-ndirs N -- number of directions in data (default is 3)\n");
  printf("\t-nphases N -- number of phases in data (default is 5)\n");
  printf("\t-2lenses -- data acquired with 2 opposing objectives\n");
  printf("\t-bessel -- data acquired with Bessel beam SIM\n");
  printf("\t-fastSIM -- data acquired with fast live SIM\n");
  printf("\t-noff N -- number of switch-off images in NL SIM data (default is 0)\n");
  printf("\t-recalcarray -- how many times do you want to re-calculuate overlapping arrays (default is 1)\n");
  printf("\t-inputapo N -- number of pixels to apodize the input data (-1 to cosine apodize)\n");
  printf("\t-forcemodamp f1 f2... f_norders -- force the modulation amplitude to be f1 and f2\n\t\tIf other than 3 phases are used, the -nphases flag must be used BEFORE the -forcemodamp flag\n");
  printf("\t-nok0search -- do not want to search for the best k0\n");
  printf("\t-nokz0 -- do not use kz0 plane in makeoverlaps() or assembly (for bad \n");
  printf("-k0searchAll -- search for k0 for every time point in a time series\n");
  printf("\t-k0angles f0 f1... f_(ndirs-1) -- user supplied pattern angle list, the -ndirs flag must be used BEFORE the -k0angles flag\n");
  printf("\t-k0spacingsf0 f1... f_(ndirs-1) -- user supplied spacings list, the -ndirs flag must be used BEFORE the -k0spacings flag\n");
  printf("\t-fitonephase -- in 2D NLSIM for orders > 1, modamp's phase will be order 1 phase multiplied by order; default is using fitted phases for all orders\n");
  printf("\t-noapodizeout -- do not want to apodize the output data\n");
  printf("\t-gammaApo -- apodize the output data with a power function\n");
  printf("\t-zoomfact -- factor by which to subdivide pixels laterally\n");
  printf("\t-zzoom -- factor by which to subdivide pixels axially\n");
  printf("\t-zpadto -- how many total sections after zero padding\n");
  printf("\t-explodefact -- factor by which to exaggerate the order shifts for display\n");
  printf("\t-nofilteroverlaps -- (Used with explodefact) leave orders round (no filtering the overlapped regions); default=do filtering \n");
  printf("\t-nosuppress -- do not want to suppress singularity at OTF origins\n");
  printf("\t-suppressR -- the radius of range\n");
  printf("\t-dampenOrder0 -- dampen order 0 in filterbands\n");
  printf("\t-noOrder0 -- do not use order 0 in assembly\n");
  printf("\t-noequalize -- no equalization of input data\n");

  printf("\t-equalizez -- to equalize images of all z sections and directions\n");
  printf("\t-equalizet -- to equalize all time points based on the first one\n");
  printf("\t-wiener f -- set wiener constant to be f (default is 0.00)\n");
  printf("\t-wienerInr f -- wiener constant of the final time point will be wiener + this number (default is 0.00)\n");
  printf("\t-background f -- set the constant background intensity (default is 515)\n");
  printf("\t-bgInExtHdr -- the background of each section is recorded in the extended header's 3rd float number (in Lin's scope)\n");
  printf("\t-otfRA -- to use radially averaged OTFs");
  printf("\t-driftfix -- to estimate and then fix drift in 3D\n");
  printf("\t-driftHPcutoff x -- the cutoff frequency (fraction of lateral resolution limit) used in high-pass Gaussian filter in determining drift (default 0.0)\n");
  printf("\t-fix2Ddrift -- to correct 2D drifts for each exposure within a pattern direction \n");
  printf("\t-fixphasestep -- to correct phase used in separation matrix based on within-direction drift correction \n");
  printf("\t-noff -- number of switch-off images in nonlinear SIM \n");
  printf("\t-usecorr corr_file -- use correction file to do flatfielding\n");
  printf("\t-nordersout x -- the number of orders to use in reconstruction. Use if different from (n_phases+1)/2\n");
  printf("\t-angle0 x -- the starting angle (in radians) of the patterns\n");
  printf("\t-negDangle -- use negative angle step\n");
  printf("\t-ls x -- the illumination pattern's line spacing (in microns)\n");
  printf("\t-na x -- the (effective) NA of the objective\n");
  printf("\t-nimm x -- the index of refraction of the immersion liquid\n");
  printf("\t-saveprefiltered filename -- save separated bands into file\n");
  printf("\t-savealignedraw filename -- save drift-corrected raw images into file\n");
  printf("\t-saveoverlaps filename -- save overlaps by makeoverlaps() into file\n");
  printf("\t-help or -h -- print this message\n");
}

int commandline(int argc, char *argv[], ReconParams *myParams, char *ifiles, char *ofiles, char * otffiles, int * ifilein, int *ofilein, int * otffilein)
{
  int ncomm=1;
  int order;

  while (ncomm < argc) {
    if (strcmp(argv[ncomm], "-angle0") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->k0startangle)) {
        printf("Invalid input for switch -angle0\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-ls") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->linespacing)) {
        printf("Invalid input %s for switch -ls\n", argv[ncomm+1]);
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-na") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->na) ) {
        printf("Invalid input for switch -na\n");
        return 0;
      }
      if ( myParams->na<0.0 || myParams->na>2.0 ) {
        printf("Invalid NA value %f. NA must lie between 0 and 2.\n", myParams->na);
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-nimm") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->nimm) ) {
        printf("Invalid input for switch -nimm\n");
        return 0;
      }
      if ( myParams->nimm<1.0 || myParams->nimm>2.0 ) {
        printf("Invalid nimm value %f. nimm must lie between 1 and 2.\n", myParams->nimm);
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-ndirs") == 0) {
      if (!getinteger(argv[ncomm+1], &myParams->ndirs)) {
        printf("Invalid input for switch -ndirs\n");
        return 0;
      }
      ncomm +=2;
    }
    else if (strcmp(argv[ncomm], "-nphases") == 0) {
      if (!getinteger(argv[ncomm+1], &myParams->nphases)) {
        printf("Invalid input for switch -nphases\n");
        return 0;
      }
      ncomm +=2;
    }
    else if (strcmp(argv[ncomm], "-phasesteps") == 0) {
      myParams->phaseSteps = malloc(3*sizeof(float));
      if (!getfloats(argv+ncomm+1, myParams->phaseSteps, 3)) {
        printf("Invalid input for switch -phasesteps\n");
        return 0;
      }
      ncomm +=4;
    }
    else if (strcmp(argv[ncomm], "-phaselist") == 0) {
      myParams->phaseList = malloc(myParams->ndirs*myParams->nphases*sizeof(float));
      if (!getfloats(argv+ncomm+1, myParams->phaseList, myParams->ndirs*myParams->nphases)) {
        printf("Invalid input for switch -phaselist\n");
        return 0;
      }
      ncomm += myParams->ndirs*myParams->nphases + 1;
    }
    else if (strcmp(argv[ncomm], "-nordersout") == 0) {
      if (!getinteger(argv[ncomm+1], &myParams->norders_output)) {
        printf("Invalid input for switch -nordersout\n");
        return 0;
      }
      ncomm +=2;
    }
    else if (strcmp(argv[ncomm], "-2lenses") == 0) {
      myParams->bTwolens = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-bessel") == 0) {
      myParams->bBessel = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-fastSIM") == 0) {
      myParams->bFastSIM = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-negDangle") == 0) {
      myParams->bNegAngleStep = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-zoomfact") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->zoomfact)) {
        printf("Invalid input for switch -zoomfact\n");
        return 0;
      }
      ncomm +=2;
    }
    else if (strcmp(argv[ncomm], "-zzoom") == 0) {
      if (!getinteger(argv[ncomm+1], &myParams->z_zoom)) {
        printf("Invalid input for switch -zzoom\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-zpadto") == 0) {
      if (!getinteger(argv[ncomm+1], &myParams->nzPadTo)) {
        printf("Invalid input for switch -zpadto\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-explodefact") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->explodefact)) {
        printf("Invalid input for switch -explodefact\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-nofilteroverlaps") == 0 ) {
      myParams->bFilteroverlaps = 0;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-recalcarray") == 0) {
      if (!getinteger(argv[ncomm+1], &myParams->recalcarrays)) {
        printf("Invalid input for switch -recalarray\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-inputapo") == 0) {
      if (!getinteger(argv[ncomm+1], &myParams->napodize)) {
        printf("Invalid input for switch -inputapo\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-forcemodamp") == 0) {
      int norders;
      if (myParams->norders_output>1)
        norders = myParams->norders_output;
      else
        norders = myParams->nphases/2+1;

      if (!getfloats(argv+ncomm+1, myParams->forceamp+1, norders-1)) {
        printf("Invalid input for switch -forcemodamp\n");
        return 0;
      }
      else
        for(order=1; order<norders; order++) {
          if (myParams->forceamp[order] <=0.0 ) {
            printf("No zero or negative value allowed for forceamp\n");
            return 0;
          }
        }
      ncomm += norders;
    }
   else if (strcmp(argv[ncomm], "-ampmin") == 0) {
      int norders;
      if (myParams->norders_output>1)
        norders = myParams->norders_output;
      else
        norders = myParams->nphases/2+1;

      if (!getfloats(argv+ncomm+1, myParams->ampmin+1, norders-1)) {
        printf("Invalid input for switch -ampmin\n");
        return 0;
      }
      ncomm += norders;
    }
   else if (strcmp(argv[ncomm], "-ampmax") == 0) {
      int norders;
      if (myParams->norders_output>1)
        norders = myParams->norders_output;
      else
        norders = myParams->nphases/2+1;

      if (!getfloats(argv+ncomm+1, myParams->ampmax+1, norders-1)) {
        printf("Invalid input for switch -ampmax\n");
        return 0;
      }
      ncomm += norders;
    }
    else if (strcmp(argv[ncomm], "-k0angles") ==0 ) {
      myParams->k0angles = malloc(myParams->ndirs * sizeof(float));
      if (!getfloats(argv+ncomm+1, myParams->k0angles, myParams->ndirs)) {
        printf("Invalid input for switch -k0angles\n");
        return 0;
      }
      ncomm += myParams->ndirs+1;
    }
    else if (strcmp(argv[ncomm], "-k0spacings") ==0 ) {
      myParams->k0spacings = malloc(myParams->ndirs * sizeof(float));
      if (!getfloats(argv+ncomm+1, myParams->k0spacings, myParams->ndirs)) {
        printf("Invalid input for switch -k0spacings\n");
        return 0;
      }
      ncomm += myParams->ndirs+1;
    }
    else if (strcmp(argv[ncomm], "-nok0search") ==0 ) {
      myParams->bSearchforvector = 0;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-k0searchAll") ==0 ) {
      myParams->bUseTime0k0 = 0;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-noapodizeout") == 0 ){
      myParams->apodizeoutput = 0;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-gammaApo") == 0) {
      myParams->apodizeoutput = 2;
      if (!getfloat(argv[ncomm+1], &myParams->apoGamma)) {
        printf("Invalid input %s for switch -gammaApo\n", argv[ncomm+1]);
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-preciseapo") == 0) {
      myParams->bPreciseApoBoundary = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-nosuppress") == 0) {
      myParams->bSuppress_singularities = 0;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-suppressR") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->suppression_radius)) {
        printf("Invalid input for switch -supressR\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-dampenOrder0") == 0) {
      myParams->bDampenOrder0 = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-noOrder0") == 0) {
      myParams->bNoOrder0 = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-fitonephase") == 0) {
      myParams->bFitallphases = 0;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-noequalize") == 0) {
      myParams->do_rescale = 0;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-otfRA") == 0) {
      myParams->bRadAvgOTF = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-equalizez") == 0) {
      myParams->equalizez = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-equalizet") == 0) {
      myParams->equalizet = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-nokz0") == 0) {
      myParams->bNoKz0 = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-2D") == 0) {
      myParams->b2D = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-wiener") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->wiener) || myParams->wiener<0 ) {
        printf("Invalid input for switch -wiener\n");
        return 0;
      }
      ncomm += 2;
      if (myParams->wiener > 0)
        myParams->bUseEstimatedWiener = 0;
    }
    else if (strcmp(argv[ncomm], "-wienerInr") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->wienerInr)) {
        printf("Invalid input for switch -wienerInr\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-driftfix") == 0) {
      myParams->bFixdrift = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-driftHPcutoff") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->drift_filter_fact)) {
        printf("Invalid input for switch -driftHPcutoff\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-fix2Ddrift") == 0) {
      myParams->bFix2Ddrift = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-fixphasestep") == 0) {
      myParams->bFixPhaseStep = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-noff") == 0) {
      if (!getinteger(argv[ncomm+1], &myParams->nOffImages)) {
        printf("Invalid input for switch -noff\n");
        return 0;
      }
      ncomm +=2;
    }
    else if (strcmp(argv[ncomm], "-fitdrift") == 0) {
      myParams->bFitDrift = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-combineoff") == 0) {
      int norders;
      if (myParams->norders_output > 1)
        norders = myParams->norders_output;
      else
        norders = myParams->nphases/2+1;

      if (!getints(argv+ncomm+1, myParams->combineOffOrders, norders)) {
        printf("Invalid input for switch -combineoff\n");
        return 0;
      }
      else
        for(order=0; order<norders; order++) {
          if (myParams->combineOffOrders[order] > 0 ) {
            myParams->bCombineOffExp = 1;
            break;
          }
        }
      ncomm += norders+1;

    }
    else if (strcmp(argv[ncomm], "-background") == 0) {
      if (!getfloat(argv[ncomm+1], &myParams->constbkgd)) {
        printf("Invalid input for switch -background\n");
        return 0;
      }
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-bgInExtHdr") == 0) {
      myParams->bBgInExtHdr = 1;
      ncomm += 1;
    }
    else if (strcmp(argv[ncomm], "-usecorr") == 0) {
      myParams->bUsecorr = 1;
      strcpy(myParams->corrfiles, argv[ncomm+1]);
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-makemodel") == 0) {
      myParams->bMakemodel = 1;
      ncomm ++;
    }
    else if (strcmp(argv[ncomm], "-saveprefiltered") == 0) {
      myParams->bSaveSeparated = 1;
      strcpy(myParams->fileSeparated, argv[ncomm+1]);
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-savealignedraw") == 0) {
      myParams->bSaveAlignedRaw = 1;
      strcpy(myParams->fileRawAligned, argv[ncomm+1]);
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-saveoverlaps") == 0) {
      myParams->bSaveOverlaps = 1;
      strcpy(myParams->fileOverlaps, argv[ncomm+1]);
      ncomm += 2;
    }
    else if (strcmp(argv[ncomm], "-nthreads") == 0) {
      if (!getinteger(argv[ncomm+1], &myParams->nthreads)) {
        printf("Invalid input for switch -nthreads\n");
        return 0;
      }
      ncomm += 2;
    }

/*     else if (strcmp(argv[ncomm], "-skip") == 0) { */
/*       if (!getinteger(argv[ncomm+1], skip_sec)) { */
/*         printf("Invalid input for switch -skip\n"); */
/*         return 0; */
/*       } */
/*       ncomm += 2; */
/*     } */
    else if (strcmp(argv[ncomm], "-help") == 0 || strcmp(argv[ncomm], "-h") == 0) {
      Usage(argv);
      return 0;
    }
    else if (ncomm < 4) {
      if (ncomm==1) {
        /* Should check whether length of file names are longer than allocated space */
        strcpy(ifiles, argv[ncomm]);
        *ifilein = 1;
      }
      else if (ncomm==2) {
        strcpy(ofiles, argv[ncomm]);
        *ofilein = 1;
      }
      else if (ncomm==3) {
        strcpy(otffiles, argv[ncomm]);
        *otffilein = 1;
      }
      ncomm ++;
    }
    else {
      printf("Invalid command line option %s\n", argv[ncomm]);
      Usage(argv);
      return 0;
    }
  }
  return 1;
}


void SetDefaultParams(ReconParams *pParams)
{
  pParams->k0startangle =1.57193;
  pParams->bNegAngleStep = 0;
  pParams->linespacing = 0.177;  /* default to Nikon TIRF 100x */
  pParams->na=1.36;
  pParams->nimm=1.515;
  pParams->ndirs = 3;
  pParams->nphases = 3;
  pParams->phaseSteps = 0;
  pParams->phaseList = 0;
  pParams->norders_output = 0;
  pParams->bTwolens = 0;
  pParams->bFastSIM = 0;
  pParams->bBessel = 0;

  pParams->zoomfact = 4;
  pParams->z_zoom = 1;
  pParams->nzPadTo = 0;
  pParams->explodefact=1.0;
  pParams->bFilteroverlaps=1;
  pParams->recalcarrays = 1; /* whether to calculate the overlaping regions between bands just once or always; used in fitk0andmodamps() */
  pParams->napodize = 10;
  pParams->forceamp[1] = 0.0;
  pParams->ampmin[1]=0.0;
  pParams->ampmax[1]=0.0;
  pParams->k0angles = NULL;
  pParams->k0spacings = NULL;
  pParams->bSearchforvector = 1;
  pParams->bUseTime0k0 = 1;  /* default to use time 0's k0 fit for the rest in a time series data */
  pParams->apodizeoutput = 1;  /* 0-no apodize; 1-cosine apodize; 2-triangle apodize; used in filterbands() */
  pParams->bPreciseApoBoundary = 0;
  pParams->bSuppress_singularities = 1;  /* whether to dampen the OTF values near band centers; used in filterbands() */
  pParams->suppression_radius = 10;   /* if suppress_singularities is 1, the range within which suppresion is applied; used in filterbands() */
  pParams->bDampenOrder0 = 0;  /* whether to dampen order 0 contribution; used in filterbands() */
  pParams->bNoOrder0 = 0;  /* whether to dampen order 0 contribution; used in filterbands() */
  pParams->bFitallphases = 1;
  pParams->do_rescale=1; /* fading correction method: 0-no correction; 1-with correction */
  pParams->equalizez = 0;
  pParams->equalizet = 0;
  pParams->bNoKz0 = 0;
  pParams->wiener = 0.0;   /* if no >0 wiener value is given, then use the estimated wiener */
  pParams->wienerInr = 0.0;   /* Used in time series: wiener+wienerInr will be the wiener used for the final time point*/
  pParams->bUseEstimatedWiener = 1;

  pParams->bRadAvgOTF = 0;  /* default to use non-radially averaged OTFs */
  pParams->b2D = 0;

  pParams->bFixdrift = 0;
  pParams->drift_filter_fact = 0.0;
  pParams->bFix2Ddrift = 0;
  pParams->bFixPhaseStep = 0;
  pParams->nOffImages = 0;
  pParams->type_of_refImage = 0; /* default to use sum of exp and off images as reference images */
  pParams->bFitDrift = 0;
  pParams->bCombineOffExp = 0;
  memset(pParams->combineOffOrders, 0, MAXPHASES * sizeof(int));

  pParams->constbkgd = 0.0;
  pParams->bBgInExtHdr = 0;
  pParams->bUsecorr = 0;
  pParams->corrfiles[0] = '\0';  /* name of CCD correction file if usecorr is 1 */
  pParams->electrons_per_bit = 0.6528;
  pParams->readoutNoiseVar = 32.42;  // electron^2

  pParams->bMakemodel = 0;
  pParams->bSaveAlignedRaw = 0;
  pParams->fileSeparated[0] = '\0';
  pParams->bSaveSeparated = 0;
  pParams->fileRawAligned[0] = '\0';
  pParams->bSaveOverlaps = 0;
  pParams->fileOverlaps[0] = '\0';

  pParams->nthreads = -1;
}

void mrc_file_write(const float *buffer, int nx, int ny, int nz, float rlen, float zlen, int mode, int iwave, const char *files)
{
  int ostream_no=19;
  IW_MRC_HEADER header;
  int dimx, dimy, dimz, nxy, i;
  float amin=0, amax=1, amean=0.1; /* make-shift way of doing this */

  printf("Writing output file: %s\n", files);

  if (IMOpen(ostream_no, files, "new")) {
    fprintf(stderr, "File %s can not be created.\n", files);
    exit(-1);
  }

  dimx = nx;
  dimy = ny;
  dimz = nz;
  nxy = nx*ny;

  switch (mode) {
  case 0: /* normal floating-point data */
    header.mode = IW_FLOAT;
    break;
  case 1: /* half Fourier space complex data */
    header.mode = IW_COMPLEX;
    dimx = nx/2 + 1;
    nxy = (nx+2)*ny;
    break;
  case 2: /* full Fourier space complex data */
    header.mode = IW_COMPLEX;
    nxy = nx*ny*2;
    break;
  case 3: /* half Fourier floating-point data */
    header.mode = IW_FLOAT;
    dimx = nx/2+1;
    nxy = dimx * dimy;
    break;
  default:
    fprintf(stderr, "Illegal mode in mrc_file_write()\n");
    exit(-1);
  }

  header.nx = dimx;
  header.mx = dimx;
  header.ny = dimy;
  header.my = dimy;
  header.nz = dimz;
  header.mz = dimz;
  header.ylen = rlen;
  header.xlen = rlen;
  header.zlen = zlen;
  header.nxst = 0;
  header.nyst = 0;
  header.nzst = 0;
  header.num_waves = 1;
  header.num_times = 1;
  header.iwav1 = iwave;
  header.alpha = 0;
  header.beta = 0;
  header.gamma = 0;
  header.mapc = 1;
  header.mapr = 2;
  header.maps = 3;
  header.ispg = 0;
  header.nDVID = -16244;
  header.ntst = 0;
  header.inbsym = 0;
  header.nint = 0;
  header.nreal = 0;
  header.nres = 1;
  header.nzfact = 1;
  header.file_type = 0;  /* normal image type */
  header.lens = 0;
  header.interleaved = 2;
  header.nlab = 1;
  IMPutHdr(ostream_no, &header);

  for (i=0; i<nz; i++)
    IMWrSec(ostream_no, buffer+i*nxy);

  IMWrHdr(ostream_no, "Processed", 1, amin, amax, amean);
  IMClose(ostream_no);
}

void load_and_flatfield(int istream_no, int section_no, int wave_no, int time_no, float *bufDestiny, float *buffer, int nx, int ny, float *background, float backgroundExtra, float *slope, float inscale, int bUsecorr)
/*
  Load the next 2D section from the MRC file identified by "istream_no".
  "bufDestiny" is where the current loaded section ends up being; it's assumed to have 2 extra columns for in-place FFT later
  "buffer" is a nx*ny sized array to hold temporarily the loaded data before it is flat-fielded and copied to "bufDestiny"
*/
{
  int l, k;

  IMPosnZWT(istream_no, section_no, wave_no, time_no);
  IMRdSec(istream_no, buffer);

  /* fix camera error at even binnings (on old OM2 CCD) */
  if(buffer[(ny-1)*nx + (nx-1)] == -1.0) {
    buffer[(ny-1)*nx + (nx-1)] = buffer[(ny-1)*nx + (nx-2)];
    slope[(ny-1)*nx + (nx-1)] = 1;
    background[(ny-1)*nx + (nx-1)] = background[(ny-1)*nx + (nx-2)];
  }

  for (l=0; l<ny; l++) {
	  for (k=0; k<nx; k++) {
		  if (bUsecorr)
			  bufDestiny[l*(nx+2) + k] = ((buffer[l*nx + k]-background[l*nx+k]-backgroundExtra) * slope[l*nx+k]) * inscale;
		  else
			  bufDestiny[l*(nx+2) + k] = (buffer[l*nx + k]-background[l*nx+k]-backgroundExtra) * inscale;
	  }
    for(k=nx;k<nx+2;k++)
      bufDestiny[l*(nx+2) + k] = 0.0;
  }

}

#ifdef __SIRECON_USE_TIFF__
int loadOTFs(CImg & otf_tiff, fftwf_complex **otf, int norders, int nz, int sizeOTF, ReconParams *pParams) {
  uint32 nxotf, nyotf, nzotf;
  float xres, yres;
  nxotf = otf_tiff.width();
  nyotf = otf_tiff.height();
  /* TIFFGetField(otf_tiff, TIFFTAG_XRESOLUTION, &xres); */
  /* TIFFGetField(otf_tiff, TIFFTAG_YRESOLUTION, &yres); */
  nzotf = otf_tiff.depth();
  if (nz == 1) {  /* 2D */
    pParams->nxotf = nxotf;
    if (pParams->bRadAvgOTF)
      pParams->nyotf = 1;
    else
      pParams->nyotf = nyotf;
    pParams->nzotf = 1;
    pParams->dkrotf = xres;  // dkrotf's unit is 1/micron
    pParams->dkzotf = 1;
  }
  else {   /* 3D */
    if (pParams->bRadAvgOTF) {
      pParams->nzotf = nxotf;
      pParams->nxotf = nyotf;
      pParams->nyotf = 1;
      pParams->dkzotf = xres;
      pParams->dkrotf = yres;
    }
    else {
      pParams->nzotf = nzotf / norders; // each order has a 3D OTF stack (non-negative kx half of Fourier space)
      pParams->nxotf = nxotf;
      pParams->nyotf = nyotf;
      pParams->dkzotf = yres; // use yres for dkzotf
      pParams->dkrotf = xres; // use xres for dkrotf
    }
  }
  printf("nzotf=%d, dkzotf=%f, nxotf=%d, nyotf=%d, dkrotf=%f\n",pParams->nzotf, pParams->dkzotf, pParams->nxotf, pParams->nyotf, pParams->dkrotf);

  /* sizeOTF and norders are determined so that correct memory can be allocated for otf */
  sizeOTF = pParams->nzotf*pParams->nxotf*pParams->nyotf;

  /* Load OTF data, no matter 2D, 3D, radially averaged or not. */
  for (i=0; i<norders; i++) {
    /* If OTF file has multiple sections, then read them into otf[i]; */
    if (nz == 1 || pParams->bRadAvgOTF) {
      if (nzotf > i)  /* each section in OTF file is OTF of one order; so load that section into otf[i]  */
        memcpy(otf[i], otf_tiff.data(0, 0, i, 0), sizeOTF * sizeof(fftwf_complex));
      else   /* If there's just 1 OTF image, do not read any more and just duplicate otf[0] into otf[i] */
        memcpy(otf[i], otf[0], sizeOTF * sizeof(fftwf_complex));
    }
    else {  // non-radially averaged 3D OTF
      memcpy(otf[i], sizeOTF * sizeof(fftwf_complex));
    }
  }
  TIFFclose(otf_tiff);
  return 1;
}
#else
void determine_otf_dimensions(int otfstream_no, int norders, int nz, ReconParams *pParams, int *sizeOTF)
{
  int ixyz[3], mxyz[3], pixeltype;
  float min, max, mean;
  IW_MRC_HEADER otfheader;
  /* Retrieve OTF file header info */
  IMRdHdr(otfstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(otfstream_no, &otfheader);
  IMAlCon(otfstream_no, 0);
  /* determine nzotf, nxotf, nyotf, dkrotf, dkzotf based on dataset being 2D/3D and flag bRadAvgOTF */
  if (nz == 1) {  /* 2D */
    pParams->nxotf = otfheader.nx;
    if (pParams->bRadAvgOTF)
      pParams->nyotf = 1;
    else
      pParams->nyotf = otfheader.ny;
    pParams->nzotf = 1;
    pParams->dkrotf = otfheader.xlen;  // dkrotf's unit is 1/micron
    pParams->dkzotf = 1;
  }
  else {   /* 3D */
    if (pParams->bRadAvgOTF) {
      pParams->nzotf = otfheader.nx;
      pParams->nxotf = otfheader.ny;
      pParams->nyotf = 1;
      pParams->dkzotf = otfheader.xlen;
      pParams->dkrotf = otfheader.ylen;
    }
    else {
      pParams->nzotf = otfheader.nz / norders; // each order has a 3D OTF stack (non-negative kx half of Fourier space)
      pParams->nxotf = otfheader.nx;
      pParams->nyotf = otfheader.ny;
      pParams->dkzotf = otfheader.zlen;
      pParams->dkrotf = otfheader.xlen;
    }
  }
  printf("nzotf=%d, dkzotf=%f, nxotf=%d, nyotf=%d, dkrotf=%f\n",pParams->nzotf, pParams->dkzotf, pParams->nxotf, pParams->nyotf, pParams->dkrotf);

  /* sizeOTF and norders are determined so that correct memory can be allocated for otf */
  *sizeOTF = pParams->nzotf*pParams->nxotf*pParams->nyotf;
}

int loadOTFs(int otfstream_no, fftwf_complex **otf, int norders, int nz, int sizeOTF, ReconParams *pParams) {
  int i, z;
  int ixyz[3], mxyz[3], pixeltype;
  float min, max, mean;
  IW_MRC_HEADER otfheader;
  /* Retrieve OTF file header info */
  IMRdHdr(otfstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(otfstream_no, &otfheader);

  /* Load OTF data, no matter 2D, 3D, radially averaged or not. */
  for (i=0; i<norders; i++) {
    /* If OTF file has multiple sections, then read them into otf[i]; */
    if (nz == 1 || pParams->bRadAvgOTF) {
      if (otfheader.nz > i)  /* each section in OTF file is OTF of one order; so load that section into otf[i]  */
        IMRdSec(otfstream_no, otf[i]);
      else   /* If there's just 1 OTF image, do not read any more and just duplicate otf[0] into otf[i] */
        memcpy(otf[i], otf[0], sizeOTF * sizeof(fftwf_complex));
    }
    else {  // non-radially averaged 3D OTF
      for (z=0; z<pParams->nzotf; z++)
        IMRdSec(otfstream_no, otf[i]+z*pParams->nxotf*pParams->nyotf);
    }
  }
  IMClose(otfstream_no);
  return 1;
}
#endif


void getbg_and_slope(char *corrfiles, float *background, float *slope, int nx, int ny)
{
  int cstream_no=10;
  int ixyz[3], mxyz[3], pixeltype;      /* variables for IMRdHdr call */
  float min, max, mean;      /* variables for IMRdHdr call */
  IW_MRC_HEADER header;

  if (IMOpen(cstream_no, corrfiles, "ro")) {
    fprintf(stderr, "File %s does not exist\n", corrfiles);
    exit(-1);
  }

  IMRdHdr(cstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(cstream_no, &header);
  if (header.nx != nx || header.ny != ny) {
    fprintf(stderr, "calibration file %s has different dimension than data file", corrfiles);
    exit(-1);
  }
  IMAlCon(cstream_no, 0);

  IMRdSec(cstream_no, background);
  IMRdSec(cstream_no, slope);

  IMClose(cstream_no);
}
