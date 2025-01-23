#include <IMInclude.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> // memset()

#include "sirecon.h"


void determinedrift_bt_dirs(int stream_no, int wave_no, int time_no, vector3d* drift, int nx, int ny, int nz, int ndirs, int nphases, float *background, float dxy, float dz, float rdistcutoff, short wave, ReconParams *pParams)
{
/* 
   Checks for slow drift between the different directions by cross-correlating 
   the (nominally identical) center bands.
*/
  float ** bandcenter, *real_crosscorr, maxval, reval, valplus, valminus;
  int i, j, k, dir, phase, nxy, nxy2, nxyz, maxi, maxj, maxk, kplus, kminus, iplus, iminus, jplus, jminus, ind;
  float * buffer;
  fftwf_complex **Cbandcenter, *crosscorr;
  fftwf_plan fftplanf, fftplanb;
  float dkx, dky, dkz, zdistcutoff;
  int subz=nz/2;    /* the number of middle sections to include in drift estimation */
  float power = 0.25;  /* power function of the crosscorr's magnitude is used */
  int nDirPhases = ndirs * nphases;

  printf("In determinedrift_bt_dirs\n\n");
  if(subz==0) subz=1;  /* In 2D case, nz/2 = 0 */
  dkx = 1/(nx*dxy);   /* inverse microns per pixel in data */
  dky = 1/(ny*dxy);   /* inverse microns per pixel in data */
  dkz = 1/(nz*dz);   /* inverse microns per pixel in data */
  /* rdistcutoff = (int)( (pParams->na*2/(wave/1000.0)) / dkr );    /\* OTF support radial limit in data pixels *\/ */
  /* if( rdistcutoff>nx/2 ) rdistcutoff=nx/2; */
  zdistcutoff = (int)ceil( ( (1-cos(asin(pParams->na/pParams->nimm)))/((wave/pParams->nimm)/1000.0) ) / dkz );    /* OTF support axial limit in data pixels */
  if( zdistcutoff>nz/2 ) zdistcutoff=nz/2;

  real_crosscorr = (float *) malloc( ((nx+2)*ny*subz) * sizeof(float));
  crosscorr = (fftwf_complex *)real_crosscorr;

  nxy = nx*ny;
  nxy2= (nx+2)*ny;
  nxyz = nz*nxy;

  bandcenter = (float **) malloc(ndirs * sizeof(float *));
  Cbandcenter = (fftwf_complex **) malloc(ndirs * sizeof(fftwf_complex *));
  for (i=0; i<ndirs; i++) {
    bandcenter[i] = (float *) malloc(nxy2*subz*sizeof(float));
    if (!bandcenter[i]) {
      printf("Not enough memory in \"determinedrift_bt_dirs\"\n");
      exit(-1);
    }
    memset(bandcenter[i], 0, nxy2*subz*sizeof(float));
    Cbandcenter[i] = (fftwf_complex *) bandcenter[i];
  }
  buffer = (float *) malloc(nxy*sizeof(float));

  if (!real_crosscorr ||!buffer || !bandcenter[0]) {
    printf("Not enough  memory in \"determinedrift_3D\"\nProgram aborted\n");
    exit(-1);
  }

  if (subz>1) {
    fftplanf = fftwf_plan_dft_r2c_3d(subz, ny, nx, bandcenter[0], Cbandcenter[0], FFTW_ESTIMATE);
    fftplanb = fftwf_plan_dft_c2r_3d(subz, ny, nx, crosscorr, real_crosscorr, FFTW_ESTIMATE);
  }
  else {
    fftplanf = fftwf_plan_dft_r2c_2d(ny, nx, bandcenter[0], Cbandcenter[0], FFTW_ESTIMATE);
    fftplanb = fftwf_plan_dft_c2r_2d(ny, nx, crosscorr, real_crosscorr, FFTW_ESTIMATE);
  }

  /* start at the middle subz sections of the first stack */
  /* current file pointer is right at the start of the first stack */
  //  IMPosn(stream_no, (nz-subz)/2*nphases*(1+nOffImages), 0);

  for (dir=0; dir<ndirs; dir++) {
    for (k=0; k<subz; k++) {
      int zsec;
      /* First determine which section of the raw data to load */
      if (pParams->bFastSIM)  /* data organized into (nz, ndirs, nphases) */
        zsec = ((k+(nz-subz)/2)*nDirPhases + dir*nphases) * (1+pParams->nOffImages);
      else   /* data organized into (ndirs, nz, nphases) */
        zsec = (dir*nz*nphases + (k+(nz-subz)/2)*nphases) * (1+pParams->nOffImages);

      for (phase=0; phase<nphases; phase++) {
        /* Load switch-off images if any is available */
        if (pParams->nOffImages>0)
          for (i=0; i<pParams->nOffImages; i++) {
            IMPosnZWT(stream_no, zsec, wave_no, time_no);
            IMRdSec(stream_no, buffer);
            if(buffer[(ny-1)*nx + (nx-1)] == -1.0)   /* fix error that occurs on OM2 Photometrics camera at even binnings */
              buffer[(ny-1)*nx + (nx-1)] = buffer[(ny-1)*nx + (nx-2)];
            for (i=0; i<ny; i++)
              for (j=0; j<nx; j++)
                bandcenter[dir][k*nxy2+i*(nx+2)+j] += (buffer[i*nx+j]-background[i*nx+j])/(float) nphases / nxyz;
          }

        zsec += pParams->nOffImages;
        

        IMPosnZWT(stream_no, zsec, wave_no, time_no);
        IMRdSec(stream_no, buffer);
        if(buffer[(ny-1)*nx + (nx-1)] == -1.0)   /* fix error that occurs on OM2 Photometrics camera at even binnings */
          buffer[(ny-1)*nx + (nx-1)] = buffer[(ny-1)*nx + (nx-2)];
        for (i=0; i<ny; i++)
          for (j=0; j<nx; j++)
            bandcenter[dir][k*nxy2+i*(nx+2)+j] += (buffer[i*nx+j]-background[i*nx+j])/(float)nphases / nxyz;
        zsec ++;
      }

      apodize(10, nx, ny, &(bandcenter[dir][k*nxy2]));
    }
    //    if (dir<ndirs-1) IMPosn(stream_no, ((nz-subz)/2+(dir+1)*nz)*nphases*(1+nOffImages), 0); /* jump to the middle of the next stack */
  }

  for (dir=0; dir<ndirs; dir++)   /* into Fourier space */
    fftwf_execute_dft_r2c(fftplanf, bandcenter[dir], Cbandcenter[dir]);

  /* Cbandcenter[] now has half Fourier space info */

  for (dir=1; dir<ndirs; dir++) {
    memset(real_crosscorr, 0, ((nxy2*subz) * sizeof(float)));
    for (k=0; k<subz; k++)
      for (i=0; i<ny; i++)
        for (j=0; j<nx/2+1; j++) {
          fftwf_complex val0, val;
          int ind, ishift;
          float rdist;

          ind = k*(nx/2+1)*ny + i*(nx/2+1) + j;
          val = Cbandcenter[dir][ind];
/*           val.im *= -1; */

          /* want to test if high-pass filtering before cross-correlation might help */
/*           kshift = k; */
          ishift = i;
/*           if (k>subz/2) kshift -= nz; */
          if (i>ny/2) ishift -= ny;
          rdist=sqrt(ishift*dky*ishift*dky + j*dkx*j*dkx);
          if (rdist<rdistcutoff) {
            val0 = conjf(Cbandcenter[0][ind]);
            crosscorr[ind] = val0 * val;
            float mag = cabsf(crosscorr[ind]), phase=cargf(crosscorr[ind]);
            crosscorr[ind] = pow(mag, power) * (cos(phase) + I*sin(phase));
            // high-pass filtering
            if (pParams->drift_filter_fact > 0.0) {
              float rho=rdist/(pParams->drift_filter_fact*rdistcutoff);
              float fact = 1-exp(-rho*rho);
              crosscorr[ind] *= fact*fact;
            }
          }
        }

    fftwf_execute(fftplanb);   /* back to real space */

    /*now crosscorr is 3d cross-correlation*/
    maxval=0.0; maxj=nx/2; maxi=ny/2; maxk=subz/2;
    for(k=0;k<subz;k++)
      for(i=0;i<ny;i++)
        for(j=0;j<nx;j++)
          {
            ind=k*nxy2+i*(nx+2)+j;
            reval=fabs(real_crosscorr[ind]);
            if( reval > maxval ) {
              maxval = reval;
              maxi=i; maxj=j;
              maxk=k;
            }
          }
    /*  locate max pixel to subpixel accuracy by fitting parabolas  */
    /*    Or, back in fourier space, do a weighted plane fit to the relative phase?  */
    iminus = maxi-1; iplus = maxi+1;
    if( iminus<0 ) iminus+=ny;
    if( iplus>=ny ) iplus-=ny;
    jminus = maxj-1; jplus = maxj+1;
    if( jminus<0 ) jminus+=nx;
    if( jplus>=nx ) jplus-=nx;
    kminus = maxk-1; kplus = maxk+1;
    if( kminus<0 ) kminus+=subz;
    if( kplus>=subz ) kplus-=subz;

    /* BE CAREFUL -- real_crosscorr is of dimension (nx+2, ny, nz) */
    valminus = fabs(real_crosscorr[ kminus*nxy2+maxi*(nx+2)+maxj ]);
    valplus  = fabs(real_crosscorr[ kplus *nxy2+maxi*(nx+2)+maxj ]);
    if(subz>1)
      drift[dir].z = maxk + fitparabola( valminus, maxval, valplus );
    else
      drift[dir].z = 0;
    if( drift[dir].z>subz/2 ) drift[dir].z-= subz;

    valminus = fabs(real_crosscorr[ maxk*nxy2+iminus*(nx+2)+maxj ]);
    valplus  = fabs(real_crosscorr[ maxk*nxy2+iplus *(nx+2)+maxj ]);
    drift[dir].y = maxi + fitparabola( valminus, maxval, valplus );
    if( drift[dir].y>ny/2 ) drift[dir].y-= ny;

    valminus = fabs(real_crosscorr[ maxk*nxy2+maxi*(nx+2)+jminus ]);
    valplus  = fabs(real_crosscorr[ maxk*nxy2+maxi*(nx+2)+jplus  ]);
    drift[dir].x = maxj + fitparabola( valminus, maxval, valplus );
    if( drift[dir].x>nx/2 ) drift[dir].x-= nx;

  }

  for (i=0; i<ndirs; i++)
    free(bandcenter[i]);
  free(bandcenter);
  free(Cbandcenter);
  free(buffer);
  free(real_crosscorr);
  fftwf_destroy_plan(fftplanf);
  fftwf_destroy_plan(fftplanb);
}



void fixdrift_bt_dirs(fftwf_complex *bands[], int norders, vector3d drift, int nx,int ny, int nz)
/*
  Applies a correcting phase factor in Fourier space, to both center and side bands;
  to correct the overall between-directions drift for this particular direction.
  Input "bands" is already in Fourier space   
*/
{
  fftwf_complex outval, expiphi;
  int order, i, j, k, ycent, zcent, nxy, kz, ky, kx;
  float angle;

  printf("In fixdrift_bt_dirs()\n\n");
  /* apply the phase correction factor */
  ycent = ny/2; zcent=nz/2;
  nxy = (nx/2+1)*ny;

  for (order=0; order<norders; order++) {
#pragma omp parallel for shared(bands, order, ycent, zcent, nxy) private(k, i, j, kz, kx, ky, angle, expiphi, outval)
    for(k=0;k<nz;k++) {
      int indz=k*nxy;
      kz=k; if(k>zcent) kz-=nz;
      for(i=0;i<ny;i++) {
        int indzy = indz + i*(nx/2+1);
        ky=i; if(i>ycent) ky-=ny;
        for(j=0;j<nx/2+1;j++) {
          int ind;
          kx=j;
          angle = 2* M_PI* ( drift.x*kx/nx + drift.y*ky/ny + drift.z*kz/nz );
          expiphi = cos(angle) + I*sin(angle);
          ind = indzy + j;
          if(order==0) {
            outval = bands[0][ind] * expiphi;
            bands[0][ind] = outval; }
          else {
            outval = bands[2*order-1][ind] * expiphi;
            bands[2*order-1][ind] = outval;
            outval = bands[2*order][ind] * expiphi;
            bands[2*order][ind] = outval;
          }
        }
      }
    }
  }

}


int fitXYdrift(vector3d *drifts, float * timestamps, int nPoints, vector3d *fitted_drift, float *eval_timestamps, int nEvalPoints)
/* 
   fit a curve using the data points of (x,y) vectors in "drifts" versus "timestamps". Then evaluate drifts for 
   the nEvalPoints points using the fitted curve; return the result in fitted_drifts.
   Note: if nPoints is greater than 3, the use cubic fit? otherwise use parabola fit?
   drifts[0] should always be (0, 0)?
*/
{
  int   i, j;
  int   nPolyTerms = 4; // highest order term is nPolyTerms-1; this could be passed in as an argument
  float *matrix_A, *vecRHS;
  float *workspace, wkopt;
  int   info, lwork, nRHS=2;  // 2 because of solving for both x and y
  int   bForceZero = 1;

  if (bForceZero) {  /* for the curve to go through time point 0 */
    int nPoints_minus1, nPolyTerms_minus1;
    /* subtract all time stamps by the first */
    for (i=1; i<nPoints; i++)
      timestamps[i] -= timestamps[0];

    /* Therefore, the solution to the constant is drifts[0].
       Now we have 2 coefficients to solve. */ 

    // construct the forward matrix; it's a transpose of the C-style matrix
    matrix_A = malloc((nPoints-1) * (nPolyTerms-1) * sizeof(float));
    // leading dimension is nPoints (column in Fortran, row in C)
    for (i=0; i<nPolyTerms-1; i++)
      for (j=0; j<nPoints-1; j++)
        matrix_A[i*(nPoints-1) + j] = pow(timestamps[j+1], i+1);

    // construct the right-hand-side vectors
    vecRHS = malloc((nPoints-1) * nRHS * sizeof(float));
    for (j=0; j<nPoints-1; j++) {
      vecRHS[j]           = drifts[j+1].x - drifts[0].x;
      vecRHS[j+nPoints-1] = drifts[j+1].y - drifts[0].y;
    }

    /* Because of calling Fortran subroutines, every parameter has to be passed as a pointer */
    nPoints_minus1 = nPoints - 1;
    nPolyTerms_minus1 = nPolyTerms - 1;
    /* Query and allocate the optimal workspace */
    lwork = -1;
    sgels_("No transpose", &nPoints_minus1, &nPolyTerms_minus1, &nRHS, matrix_A, &nPoints_minus1,
          vecRHS, &nPoints_minus1, &wkopt, &lwork, &info );
    lwork = wkopt;
    workspace = malloc( lwork* sizeof(float) );

    sgels_("No transpose", &nPoints_minus1, &nPolyTerms_minus1, &nRHS, matrix_A, &nPoints_minus1,
          vecRHS, &nPoints_minus1, workspace, &lwork, &info);

    if (info != 0)
      return 0;
  
    /* Now subtract timestamp[0] from all eval_timestamps */
    for (i=0; i<nEvalPoints; i++)
      eval_timestamps[i] -= timestamps[0];

    /* Evaluate at datapoints eval_timestamps */
    for (i=0; i<nEvalPoints; i++) {
      fitted_drift[i].x = drifts[0].x;
      fitted_drift[i].y = drifts[0].y;
      for (j=0; j<nPolyTerms-1; j++) {
        fitted_drift[i].x += vecRHS[j] * pow(eval_timestamps[i], j+1);
        fitted_drift[i].y += vecRHS[j+nPoints-1] * pow(eval_timestamps[i], j+1);
      }
    }
  }
  else {
    // construct the forward matrix; it's a transpose of the C-style matrix
    matrix_A = malloc(nPoints * nPolyTerms * sizeof(float));
    // leading dimension is nPoints (column in Fortran, row in C)
    for (i=0; i<nPolyTerms; i++)
      for (j=0; j<nPoints; j++)
        matrix_A[i*nPoints + j] = pow(timestamps[j], i);

    // construct the right-hand-side vectors
    vecRHS = malloc(nPoints * nRHS * sizeof(float));

    for (j=0; j<nPoints; j++) {
      vecRHS[j] = drifts[j].x;
      vecRHS[j+nPoints] = drifts[j].y;
    }

    /* Query and allocate the optimal workspace */
    lwork = -1;
    sgels_("No transpose", &nPoints, &nPolyTerms, &nRHS, matrix_A, &nPoints, vecRHS, &nPoints, &wkopt, &lwork, &info );
    lwork = wkopt;
    workspace = malloc( lwork* sizeof(float) );

    sgels_("No transpose", &nPoints, &nPolyTerms, &nRHS, matrix_A, &nPoints, vecRHS, &nPoints, workspace, &lwork, &info);

    if (info != 0)
      return 0;
  
    /* Evaluate at datapoints eval_timestamps */
    for (i=0; i<nEvalPoints; i++) {
      fitted_drift[i].x = 0;
      fitted_drift[i].y = 0;
      for (j=0; j<nPolyTerms; j++) {
        fitted_drift[i].x += vecRHS[j] * pow(eval_timestamps[i], j);
        fitted_drift[i].y += vecRHS[j+nPoints] * pow(eval_timestamps[i], j);
      }
    }
  }

  free(matrix_A);
  free(vecRHS);
  free(workspace);

  return 1;
}


void translate2D(fftwf_complex *Cimage2D, int nx, int ny, vector3d shift)
/*
  Called by fixdrift_2D(); do translation via multiplying plane complex wave in Fourier space
*/
{
  int i, j, ycent = ny/2;

  /* Apply the plane complex wave */
  for(i=0; i<ny; i++) {
    int indy = i*(nx/2+1);
    int kx, ky;
    ky = i; 
    if (i>ycent)
      ky -= ny;
    for (j=0; j<nx/2+1; j++) {
      int ind;
      fftwf_complex exp_i_phi;
      float angle;
      kx = j;
      angle = 2* M_PI* ( shift.x * kx/nx + shift.y * ky/ny);
      exp_i_phi = cos(angle) + I*sin(angle);
      ind = indy + j;
      Cimage2D[ind] *= exp_i_phi;

      /* scale by 1/(nx*ny) */
/*       Cimage2D[ind].re /= nxy; */
/*       Cimage2D[ind].im /= nxy; */
    }
  }

/*   /\* Create 2D backward RFFTW plan if rfftplan2B is not provided *\/ */
/*   if (rfftplan2B == NULL) */
/*     rfftplan2B = rfftw2d_create_plan(ny, nx, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE); */

/*   /\* Backward FFT Cimage2D *\/ */
/*   rfftwnd_threads_one_complex_to_real(nthreads, rfftplan2B, Cimage2D, NULL); */
}

/*
"floatimage" holds all phases and all z sections for 1 direction raw, still mixed data; this function is working on only section z.
 "driftlist" is longer than nphases; the function locates fron "driftlist" the drift estimates based on dir and z
 */
void fixdrift_2D(fftwf_complex **CrawImages, vector3d *driftlist, int nphases, int nx, int ny, int nz, int dir, int z)
{
  int phase, offset;
  int nxy, nSecs_per_dir;
  
  nxy = (nx+2)*ny/2;

  nSecs_per_dir = nz * nphases;

  offset = dir*nSecs_per_dir+z*nphases;
  for (phase=0; phase<nphases; phase++) {
    translate2D(CrawImages[phase]+z*nxy, nx, ny, driftlist[offset+phase]);
  }
}


int calcRefImage(float **rawImages, float *refImage, float **offImages, int nOffImages, int nx, int ny, int nphases, int type_of_refImage)
/*
  Assuming nz=1; otherwise wrong
*/
{
  int i, j, ind2;
  int phase;

  switch (type_of_refImage) {
  case 0: /* sum of switch-off and normal exposure images of the first phase within each direction */
    phase = 0;
    for (i=0; i<ny; i++)
      for (j=0; j<nx+2; j++) {
        ind2 = i*(nx+2) + j;
        if (nOffImages)
          refImage[ind2] = rawImages[phase][ind2] + offImages[phase][ind2];
        else
          return 0;
      }
    break;
  case 1: /* sum of switch-off and exposure images; one per direction */
    for (i=0; i<ny; i++)
      for (j=0; j<nx+2; j++) {
        ind2 = i*(nx+2) + j;
        for (phase=0; phase<nphases; phase++) {
          refImage[ind2] += rawImages[phase][ind2] + offImages[phase][ind2];
        }
      }
    break;
  case 2:
    break;
  case 3:
    break;
  default:
    fprintf(stderr, "Type of reference image %d is invalid\n", type_of_refImage);
    return 0;
  }
  return 1;
}

void determinedrift_2D(fftwf_complex **CrawImages, fftwf_complex **CoffImages, int nOffImages, fftwf_complex * CrefImage,
                      vector3d *drifts, int nphases, int nx, int ny, int dir,
                      float rdistcutoff, float dkx, float dky, float drift_filter_fact)
/*
  Determine 2D drift between refImage and rawImages+offImages (if offImages are available) for each of the nphases;
  if offImages do not exist, then just use rawImages to compare with refImage.
  After between-image drift estimation, fit a curve to all the drift estimates (fitXYdrift)
 */
{
  int i, j, phase, nxy2;
  fftwf_plan fftplanb;
  float *real_crosscorr;
  fftwf_complex *CrawImage_of_one_phase,  *crosscorr;
  float power = 0.5;  /* power function of the crosscorr's magnitude used */

  nxy2 = (nx+2)*ny;

  
/*   CrefImage = malloc(nxy2/2*sizeof(fftwf_complex)); */
/*   CrawImage = malloc(nxy2/2*sizeof(fftwf_complex)); */
  crosscorr = fftwf_malloc(nxy2/2*sizeof(fftwf_complex));
  CrawImage_of_one_phase = malloc(nxy2/2*sizeof(fftwf_complex));
  real_crosscorr = malloc(nx*ny*sizeof(float));

  fftplanb = fftwf_plan_dft_c2r_2d(ny, nx, crosscorr, real_crosscorr, FFTW_ESTIMATE);

  for (phase=0; phase<nphases; phase++) {
    float maxval, reval, valminus, valplus;
    int maxi, maxj, iminus, iplus, jminus, jplus;

    for (i=0; i<ny; i++)
      for (j=0; j<nx/2+1; j++) {
        int ind = i*(nx/2+1)+j;
        CrawImage_of_one_phase[ind] = CrawImages[phase][ind];
        if (nOffImages)
            CrawImage_of_one_phase[ind] += CoffImages[phase][ind];
      }

    memset(crosscorr, 0, nxy2/2 * sizeof(fftwf_complex));

    for (i=0; i<ny; i++)
      for (j=0; j<nx/2+1; j++) {
        fftwf_complex val0, val;
        int ind, ishift;
        float rdist;

        ind = i*(nx/2+1) + j;
        val = CrawImage_of_one_phase[ind];

        ishift = i;
        if (i>ny/2) ishift -= ny;

        rdist=sqrt(ishift*dky*ishift*dky + j*dkx*j*dkx);
        if (rdist<rdistcutoff) {
          val0 = conjf(CrefImage[ind]);
          crosscorr[ind] = val0 * val;
          float mag = cabsf(crosscorr[ind]), phase=cargf(crosscorr[ind]);
          crosscorr[ind] = pow(mag, power) * (cos(phase) + I*sin(phase));
          // high-pass filtering
          if (drift_filter_fact > 0.0) {
            float rho=rdist/(drift_filter_fact*rdistcutoff);
            float fact = 1-exp(-rho*rho);
            crosscorr[ind] *= fact*fact;
          }
        }
    
      }

    fftwf_execute(fftplanb);   /* back to real space */

    maxval=0.0;
    maxj=nx/2;
    maxi=ny/2;
    for(i=0; i<ny; i++)
      for(j=0; j<nx; j++) {
        reval=fabs(real_crosscorr[i*nx + j]);
        if( reval > maxval ) {
          maxval = reval;
          maxi=i;
          maxj=j;
        }
      }
    /*  locate max pixel to subpixel accuracy by fitting parabolas  */
    /*    Or, back in fourier space, do a weighted plane fit to the relative phase?  */
    iminus = maxi-1;
    iplus = maxi+1;
    if( iminus<0 ) iminus+=ny;
    if( iplus>=ny ) iplus-=ny;
    jminus = maxj-1;
    jplus = maxj+1;
    if( jminus<0 ) jminus+=nx;
    if( jplus>=nx ) jplus-=nx;

    valminus = fabs(real_crosscorr[iminus*nx+maxj]);
    valplus  = fabs(real_crosscorr[iplus *nx+maxj]);
    drifts[phase].y = maxi + fitparabola(valminus, maxval, valplus);
    if( drifts[phase].y>ny/2 ) drifts[phase].y-= ny;

    valminus = fabs(real_crosscorr[maxi*nx+jminus]);
    valplus  = fabs(real_crosscorr[maxi*nx+jplus]);
    drifts[phase].x = maxj + fitparabola(valminus, maxval, valplus);
    if( drifts[phase].x>nx/2 ) drifts[phase].x-= nx;
  }

/*   free(CrefImage); */
  free(CrawImage_of_one_phase);
  free(crosscorr);
  free(real_crosscorr);
/*   fftwnd_destroy_plan(fftplanf); */
  fftwf_destroy_plan(fftplanb);
}

void calcPhaseList(float * phaseList, vector3d *driftlist, float *phaseAbs, float k0angle, float linespacing, float dr, int nOffImages, int nphases, int nz, int direction, int z)
/*
  Calculate the actual pattern phase value for each exposure, taking into account drift estimation, acquisition-time
  phase values, and the k0 vector of this direction.
  phaseAbs -- acquisition-time absolute phase used, recorded in extended header
*/
{
  int p;
  float phistep, k0mag;
  vector3d kvec, drift;

  phistep = 2*M_PI / nphases;

  k0mag = 1/linespacing;
  kvec.x = k0mag * cos(k0angle);
  kvec.y = k0mag * sin(k0angle);

  /* The "standard" evenly-spaced phases + phase correction caused by drift correction - 
   acquisition-time phase correction */
  for (p=0; p<nphases; p++) {
    float phiDrift;
    drift = driftlist[direction*nz*nphases + z*nphases + p];
    phiDrift = (drift.x*kvec.x + drift.y*kvec.y) * dr * 2 * M_PI;

    if (phaseAbs && nz == 1)
      // here the absolute phases of the off images are used; not necessary but simply because phase starts at 0 for off images
      phaseList[p] = phaseAbs[direction*nphases*(1+nOffImages) + p*(1+nOffImages)]  + phiDrift;
    
    else
      phaseList[p] = p * phistep + phiDrift;

    printf("phaseAbs[%d]=%10.5f, phaseList[%d]=%10.5f, phase error=%8.3f, phase drift=%8.3f deg\n",
           p, phaseAbs[direction*nphases*(1+nOffImages) + p*(1+nOffImages)]*180/M_PI,
           p, phaseList[p]*180/M_PI, 
           (phaseList[p] - p *phistep)*180/M_PI,
           phiDrift*180/M_PI);
  }

}
