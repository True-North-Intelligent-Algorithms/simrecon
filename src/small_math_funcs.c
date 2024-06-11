#include <math.h>
#include <string.h>
#include <omp.h>   // OpenMP library for multithreading via compiler
#include "sirecon.h"
#include <stdlib.h>

float mag2(fftwf_complex b)
{
  float re, im;
  re=crealf(b);
  im = cimagf(b);
  return re*re+im*im;
}


int image_arithmetic(float * imDest, float *imSrc, int len, float fDest, float fSrc)
/*
  returns fDest*imDest + fSrc*imSrc in imDest
  imSrc can be a NULL pointer, in which case only fDest*imDest will be returned
*/
{
  int i;
  if (!imDest)
    return 0;

#pragma omp parallel for private(i)
  for (i=0; i<len; i++) {
    imDest[i] *= fDest;
    if (imSrc)
      imDest[i] += fSrc * imSrc[i];
  }

  return 1;
}

double sum_of_pixels(float *im, int nx, int ny, int bExtra2Cols)
{
  int i, j, xdim;
  double sum=0;
  float *pIm;

  xdim = nx;
  if (bExtra2Cols)
    xdim += 2;

  pIm = im;
  for (i=0; i<ny; i++) {
    for (j=0; j<nx; j++)
      sum += pIm[j];
    pIm += xdim;
  }
  return sum;
}


double sqrsum_of_pixels(float *im, int nx, int ny, int bExtra2Cols)
{
  int i, j, xdim;
  double sum=0;
  float *pIm;

  xdim = nx;
  if (bExtra2Cols)
    xdim += 2;

  pIm = im;
  for (i=0; i<ny; i++) {
    for (j=0; j<nx; j++)
      sum += pIm[j]*pIm[j];
    pIm += xdim;
  }

  return sum;
}


/***************************** separate **************************************/
/*     Applies image arithmetic to the image sequence rawImages[][],        */
/*     which was acquired with different phases of the illumination          */
/*     pattern, to separate the different bands of sample information.       */
/*     The coefficients are pre-stored in the matrix sepMatrix.              */
/*     The bands are returned in the same array rawImages where the         */
/*     input data was.                                                       */
/*****************************************************************************/
void separate(int nx, int ny, int z, int direction, int nphases, int norders, float *rawImages[], float *sepMatrix)
{
  int i, j, k, l, ind, nxy2;
  float *output;

  output = (float *) malloc((norders*2-1) * sizeof(float));
  nxy2 = (nx+2)*ny;

  ind = z * nxy2;
  for (l=0; l<ny; l++) {
    for (k=0; k<nx; k++) {
      for (i=0; i<norders*2-1; i++) {
        output[i] = 0.0;
        for (j=0;j<nphases;j++) {
          output[i] += rawImages[j][ind] * sepMatrix[i*nphases+j];
        }
      }
      for (i=0;i<norders*2-1;i++)
        rawImages[i][ind] = output[i];

      /* advance ind to the next pixel (assuming row-major contiguous storage for rawImages */
      ind ++;
    }
    ind += 2;  // to account for the 2 extra columns
  }

  free(output);
}


void combine_exp_off_F(fftwf_complex **expImages, fftwf_complex **offImages, int nx, int ny, int norders, int dir, int rdistcutoff, float dkx, float dky, float *noiseVarFactors, int *combineOffOrders)
{
  float rdist2 = rdistcutoff * rdistcutoff, kr2;
  int counter_outside=0, counter_inside=0;
  double sum_exp_sqr_outside=0.0, sum_exp_sqr_inside=0.0;
  double sum_off_sqr_outside=0.0, sum_off_sqr_inside=0.0;
  float noiseVar_exp, signal_exp, noiseVar_off, signal_off, signal2_off, signal2_exp, denominator, exp_scale, off_scale;
  static float noiseVarFactor_0_0;
  int nxy2 = (nx+2)*ny, order, i, j, kx, ky, ind;

  // estimate noise variance of both expImages and offImages from outside of OTF supoort
  for (order=0; order<norders; order++) {
    if(dir==0 && order==0)
      noiseVarFactor_0_0 = noiseVarFactors[0];
    if (combineOffOrders[order] <= 0) continue;

    counter_outside=0; counter_inside=0;
    sum_exp_sqr_outside=0.0; sum_exp_sqr_inside=0.0;
    sum_off_sqr_outside=0.0; sum_off_sqr_inside=0.0;

    ind = 0;
    for (i=0; i<ny; i++) {
      ky = i;
      if (ky>ny/2)
        ky -= ny;
      ky *= dky;
      for (j=0; j<nx/2+1; j++) {
        kx = j*dkx;
        kr2 = kx*kx + ky*ky;
        if (kr2 > rdist2) {
          counter_outside += 1;
          if (order==0) {
            sum_exp_sqr_outside += mag2(expImages[order][ind]);
            sum_off_sqr_outside += mag2(offImages[order][ind]);
          }
          else {
            sum_exp_sqr_outside += mag2(expImages[2*order][ind]);
            sum_off_sqr_outside += mag2(offImages[2*order][ind]);
            sum_exp_sqr_outside += mag2(expImages[2*order-1][ind]);
            sum_off_sqr_outside += mag2(offImages[2*order-1][ind]);
          }
        }
        else {
          counter_inside += 1;
          if (order==0) {
            sum_exp_sqr_inside += mag2(expImages[order][ind]);
            sum_off_sqr_inside += mag2(offImages[order][ind]);
          }
          else {
            sum_exp_sqr_inside += mag2(expImages[2*order][ind]);
            sum_off_sqr_inside += mag2(offImages[2*order][ind]);
            sum_exp_sqr_inside += mag2(expImages[2*order-1][ind]);
            sum_off_sqr_inside += mag2(offImages[2*order-1][ind]);
          }
        }
        ind ++;
      }
    }
    noiseVar_exp = sum_exp_sqr_outside / counter_outside;
    signal2_exp = sum_exp_sqr_inside / counter_inside - noiseVar_exp;
/*     if (signal2_exp<0) {  // TODO: what to do if this happens?? */
/*       printf("WARNING: signal2_exp is less than 0 for dir=%d, order=%d\nNo combining off image\n", dir, order); */
/*       continue; */
/*     } */
    if(signal2_exp < 0)  signal2_exp = 0;

    noiseVar_off = sum_off_sqr_outside / counter_outside;
    signal2_off = sum_off_sqr_inside / counter_inside - noiseVar_off;
    if(signal2_off < 0)  signal2_off = 0;

    if(signal2_exp == 0 && signal2_off == 0)  signal2_exp = signal2_off = 1.0;

    signal_exp = sqrt(signal2_exp);
    signal_off = sqrt(signal2_off);

/*     off_scale = (signal_off / signal_exp)*(noiseVar_exp / noiseVar_off); */
/*     denominator = 1 + (signal2_off / signal2_exp)*(noiseVar_exp / noiseVar_off); */
    off_scale = signal_off  / noiseVar_off;
    exp_scale = signal_exp / noiseVar_exp;
    denominator = exp_scale + off_scale;
    //    denominator = signal2_exp / noiseVar_exp + signal2_off / noiseVar_off;
    //    denominator = 1 + (signal2_off / signal2_exp)*(noiseVar_exp / noiseVar_off);

    printf("signal_off=%e, signal_exp=%e\n", signal_off, signal_exp);
    printf("noiseVar_off=%e, noiseVar_exp=%e\n", noiseVar_off, noiseVar_exp);
    printf("exp_scale=%f, off_scale=%f, denominator=%f\n", exp_scale, off_scale, denominator);
    printf("rel_exp_scale=%f, rel_off_scale=%f\n", exp_scale/denominator, off_scale/denominator); printf("\n");
    // Combine
/*     if (order == 0) */
/*       image_arithmetic((float*)expImages[order],     (float*)offImages[order],     nxy2, 1/denominator,    off_scale/denominator); */
/*     else { */
/*       image_arithmetic((float*)expImages[2*order-1], (float*)offImages[2*order-1], nxy2, 1/denominator, -1*off_scale/denominator); */
/*       image_arithmetic((float*)expImages[2*order],   (float*)offImages[2*order],   nxy2, 1/denominator, -1*off_scale/denominator); */
/*     } */

/*     /\* Modify the list of relative noise variance to reflect the improved signal-to-noise after combining bands *\/     */
/*     noiseVarFactors[dir*norders + order] = ( noiseVar_exp + off_scale*off_scale*noiseVar_off ) / (denominator*denominator); */
/*     if(dir==0 && order==0) */
/*       noiseVarFactor_0_0 = noiseVarFactors[0]; */
/*     noiseVarFactors[dir*norders + order] /= noiseVarFactor_0_0;  // This overall normalization affects the normalization of the Wiener constant */

    if (order == 0)
      image_arithmetic((float*)expImages[order],     (float*)offImages[order],     nxy2, exp_scale/denominator,    off_scale/denominator);
    else {
      image_arithmetic((float*)expImages[2*order-1], (float*)offImages[2*order-1], nxy2, exp_scale/denominator, -1*off_scale/denominator);
      image_arithmetic((float*)expImages[2*order],   (float*)offImages[2*order],   nxy2, exp_scale/denominator, -1*off_scale/denominator);
    }

    /* Modify the list of relative noise variance to reflect the improved signal-to-noise after combining bands */
    // following makes very little noiseVarFactors; doesn't really work; TODO
    /* noiseVarFactors[dir*norders + order] = ( exp_scale*exp_scale*noiseVar_exp + off_scale*off_scale*noiseVar_off ) / (denominator*denominator); */
    /* noiseVarFactors[dir*norders + order] /= noiseVarFactor_0_0;  // This overall normalization affects the normalization of the Wiener constant */
    /* printf ("\nnoiseVarFactor_0_0=%f\n", noiseVarFactor_0_0); */
    /* printf ("noiseVarFactors[dir*norders + order]=%f\n", noiseVarFactors[dir*norders + order]); */
  }

}

void combine_exp_off(float **expImages, float **offImages, int nx, int ny, int nphases, int dir,
                     float* noiseVarFactors, float readoutNoiseVar, float electrons_per_bit)
/*
  In non-linear SIM, combine the normal exposure and the Off image, in a way where noise contributes equally to the sum
*/
{
  double sum_exp, sum_off; // sum of all real space pixel intensities (calculate for order zero only)
  double sqrsum_exp, sqrsum_off; // sum of squares of all real space pixel intensities  (calculate for each order)
  double noisevar_exp, noisevar_off, sig_to_noise2_exp, sig_to_noise2_off, scale;
  int order, norders, nxy, nxy2;
  //ToDo

  sum_exp = sum_of_pixels(expImages[0], nx, ny, 1);
  sum_off = sum_of_pixels(expImages[0], nx, ny, 1);

  norders = nphases/2 + 1;
  nxy = nx*ny;
  nxy2 = (nx+2) * ny;


  for (order=0; order<norders; order++) {
    sqrsum_exp = sqrsum_of_pixels(expImages[2*order], nx, ny, 1);
    sqrsum_off = sqrsum_of_pixels(offImages[2*order], nx, ny, 1);
    if (order>0) {
      sqrsum_exp += sqrsum_of_pixels(expImages[2*order-1], nx, ny, 1);
      sqrsum_off += sqrsum_of_pixels(offImages[2*order-1], nx, ny, 1);
    }
    noisevar_exp = (nxy * nphases * readoutNoiseVar + sum_exp * electrons_per_bit) * noiseVarFactors[dir*norders+order];
    noisevar_off = (nxy * nphases * readoutNoiseVar + sum_off * electrons_per_bit) * noiseVarFactors[dir*norders+order];
    sig_to_noise2_exp = sqrsum_exp/noisevar_exp - 1.0;
    sig_to_noise2_off = sqrsum_off/noisevar_off - 1.0;
    scale = sig_to_noise2_off / sig_to_noise2_exp;

    // Combine
    if (order == 0)
      image_arithmetic(expImages[order], offImages[order], nxy2, 1, scale);
    else {
      image_arithmetic(expImages[2*order-1], offImages[2*order-1], nxy2, 1, -1*scale);
      image_arithmetic(expImages[2*order],   offImages[2*order],   nxy2, 1, -1*scale);
    }

    /* Modify the list of relative noise variance to reflect the improved signal-to-noise after combining bands */
    noiseVarFactors[dir*norders + order] /= 1 + scale;
  }
}


/***************************** findpeak *************************************/
/*     Locates the peak of the size*size array to sub-pixel precision       */
/*     by fitting a parabola to the highest pixel and its 4 neighbors.       */
/****************************************************************************/

void findpeak(float array[], int sizex, int sizey, vector *peak)
{
  int   xcent=0, ycent=0, i, j;
  float a1, a2, a3, big;

  big = -1e11;
  for(i=0;i<sizey;i++)
    for(j=0;j<sizex;j++)
      if(array[i*sizex+j] > big) {
        big=array[i*sizex+j];
        ycent = i;  xcent = j;
      }

  if(xcent==0)
    a1 = array[ ycent*sizex +  xcent-1+sizex];
  else
    a1 = array[ ycent*sizex +  xcent-1];
  a2 = array[ ycent*sizex +  xcent  ];
  if(xcent==sizex-1)
    a3 = array[ ycent*sizex +  xcent+1-sizex];
  else
    a3 = array[ ycent*sizex +  xcent+1];
  (*peak).x = fitparabola(a1,a2,a3) + xcent;

  if(ycent==0)
    a1 = array[ (ycent-1+sizey)*sizex + xcent ];
  else
    a1 = array[ (ycent-1)*sizex + xcent ];
  a2 = array[ (ycent  )*sizex + xcent ];
  if(ycent==sizey-1)
    a3 = array[ (ycent+1-sizey)*sizex + xcent ];
  else
    a3 = array[ (ycent+1)*sizex + xcent ];
  (*peak).y = fitparabola(a1,a2,a3) + ycent;
}


/***************************** fitparabola **********************************/
/*     Fits a parabola to the three points (-1,a1), (0,a2), and (1,a3).     */
/*     Returns the x-value of the max (or min) of the parabola.             */
/****************************************************************************/

float fitparabola( float a1, float a2, float a3 )
{
  float slope,curve,peak;

  slope = 0.5* (a3-a1);         /* the slope at (x=0).  */
  curve = (a3+a1) - 2*a2;       /* (a3-a2)-(a2-a1). The change in slope per unit of x. */
  if( curve == 0 ) {
    printf("no peak: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f\n",a1,a2,a3,slope,curve);
    return( 0.0 );
  }
  peak = -slope/curve;          /* the x value where slope = 0  */
  if( peak>1.5 || peak<-1.5 ) {
    printf("bad peak position: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f, peak=%f\n",a1,a2,a3,slope,curve,peak);
    return( 0.0 );
  }
  return( peak );
}


/***************************** fitxyparabola **********************************/
/*     Fits a parabola to the three points (x1,y1), (x2,y2), and (x3,y3).     */
/*     Returns the x-value of the max (or min) of the parabola.               */
/******************************************************************************/

float fitxyparabola( float x1, float y1, float x2, float y2, float x3, float y3 )
{
  float slope1,slope2,curve,peak,xbar1,xbar2;

  if( x1==x2 || x2==x3 || x3==x1 ) {
    printf("Fit fails; two points are equal: x1=%f, x2=%f, x3=%f\n",x1,x2,x3);
    return( 0.0 );
  }
  xbar1 = 0.5 * (x1 + x2);               /* middle of x1 and x2 */
  xbar2 = 0.5 * (x2 + x3);               /* middle of x2 and x3 */
  slope1 = (y2-y1)/(x2-x1);    /* the slope at (x=xbar1).  */
  slope2 = (y3-y2)/(x3-x2);    /* the slope at (x=xbar2).  */
  curve = (slope2-slope1) / (xbar2-xbar1);       /* The change in slope per unit of x. */
  if( curve == 0 ) {
    printf("Fit fails; no curvature: r1=(%f,%f), r2=(%f,%f), r3=(%f,%f) slope1=%f, slope2=%f, curvature=%f\n",
           x1,y1,x2,y2,x3,y3, slope1,slope2,curve);
    return( 0.0 );
  }

  peak = xbar2 - slope2/curve;          /* the x value where slope = 0  */

  return( peak );
}

/*************************** suppress ************************/
/*  An optional suppression factor for the region of the     */
/*  bands near the origin, where they vary rapidly and could */
/*  cause problems related to interpolation. Used in         */
/*  filterbands.                                             */
/*************************************************************/
float suppress(float x, float rad)
{
  float x6,out;

  x6 = x*x;
  /* x6 *= x6; */

  /* out = 1.0/(1+20000/(x6+20)); */
  out = 1.0/(1+1000/(x6+.2));

  /* out = 1-exp(-0.01*x/rad) + 0.0001; */
  return out;
}


fftwf_complex otfinterpolate(fftwf_complex * otf, float kx, float ky, int kz, float kzscale, ReconParams* pParams)
/* (kx, ky, kz) is Fourier space coords with origin at kx=ky=kz=0 and going  betwen -nx(or ny,nz)/2 and +nx(or ny,nz)/2 */
{
  fftwf_complex otfval;

  if (pParams->bRadAvgOTF) {
    int irindex, izindex, indices[2][2];
    float krindex, kzindex;
    float ar, az;

    krindex = sqrt(kx*kx+ky*ky) / pParams->dkrotf;
    kzindex = kz * kzscale;
    if (kzindex<0) kzindex += pParams->nzotf;

    irindex = floor(krindex);
    izindex = floor(kzindex);

    ar = krindex - irindex;
    az = kzindex - izindex;  // az is always 0 for 2D case, and it'll just become a 1D interp

    if (izindex == pParams->nzotf-1) {
      indices[0][0] = irindex*pParams->nzotf+izindex;
      indices[0][1] = irindex*pParams->nzotf+0;
      indices[1][0] = (irindex+1)*pParams->nzotf+izindex;
      indices[1][1] = (irindex+1)*pParams->nzotf+0;
    }
    else {
      indices[0][0] = irindex*pParams->nzotf+izindex;
      indices[0][1] = irindex*pParams->nzotf+(izindex+1);
      indices[1][0] = (irindex+1)*pParams->nzotf+izindex;
      indices[1][1] = (irindex+1)*pParams->nzotf+(izindex+1);
    }
    otfval = (1-ar)*(otf[indices[0][0]]*(1-az) + otf[indices[0][1]]*az) +
      ar*(otf[indices[1][0]]*(1-az) + otf[indices[1][1]]*az);
  }
  else { /* non-radially averaged OTF */
    int ixindex, iyindex, izindex, iyindex_plus_1, izindex_plus_1, indices[2][2][2];
    float kxindex, kyindex, kzindex;
    float ax, ay, az;
    int nxyotf;
    int conj = 0;

    kxindex = kx / pParams->dkrotf;
    kyindex = ky / pParams->dkrotf;
    kzindex = kz * kzscale;
    if (kxindex<0) {
      kxindex *= -1;
      kyindex *= -1;
      kzindex *= -1;
      conj = 1;
    }
    if (kyindex<0)
      kyindex += pParams->nyotf;
    if (kzindex<0)
      kzindex += pParams->nzotf;

    ixindex = floor(kxindex);
    iyindex = floor(kyindex);
    izindex = floor(kzindex);
    ax = kxindex - ixindex;
    ay = kyindex - iyindex;
    az = kzindex - izindex;

    iyindex_plus_1 = (iyindex+1) % pParams->nyotf;
    izindex_plus_1 = (izindex+1) % pParams->nzotf;

    nxyotf = pParams->nxotf * pParams->nyotf;

    /* Find the 8 vertices surrounding the point (kxindex, kyindex, kzindex) */
    indices[0][0][0] = izindex*nxyotf        + iyindex*pParams->nxotf        + ixindex;
    indices[0][0][1] = izindex*nxyotf        + iyindex*pParams->nxotf        + (ixindex+1);
    indices[0][1][0] = izindex*nxyotf        + iyindex_plus_1*pParams->nxotf + ixindex;
    indices[0][1][1] = izindex*nxyotf        + iyindex_plus_1*pParams->nxotf + (ixindex+1);
    indices[1][0][0] = izindex_plus_1*nxyotf + iyindex*pParams->nxotf        + ixindex;
    indices[1][0][1] = izindex_plus_1*nxyotf + iyindex*pParams->nxotf        + (ixindex+1);
    indices[1][1][0] = izindex_plus_1*nxyotf + iyindex_plus_1*pParams->nxotf + ixindex;
    indices[1][1][1] = izindex_plus_1*nxyotf + iyindex_plus_1*pParams->nxotf + (ixindex+1);

    otfval = (1-az)*(otf[indices[0][0][0]]*(1-ay)*(1-ax) + otf[indices[0][0][1]]*(1-ay)*ax +
                        otf[indices[0][1][0]]*ay*(1-ax)     + otf[indices[0][1][1]]*ay*ax)
      + az*(otf[indices[1][0][0]]*(1-ay)*(1-ax) + otf[indices[1][0][1]]*(1-ay)*ax +
            otf[indices[1][1][0]]*ay*(1-ax)     + otf[indices[1][1][1]]*ay*ax);

    if (conj)
      otfval = conjf(otfval);
  }

  return otfval;
}


void matrix_transpose(float * mat, int nRows, int nCols)
{
  int i, j;
  float *tmpmat = malloc(nRows*nCols*sizeof(float));

  for (i=0; i<nRows; i++)
    for (j=0; j<nCols; j++)
      tmpmat[j*nRows+i] = mat[i*nCols+j];
  memcpy(mat, tmpmat, nRows*nCols*sizeof(float));
  free(tmpmat);
}

void print_matrix( char* desc, int m, int n, float* a, int lda )
/* Auxiliary routine: printing a matrix */
{
  int i, j;
  printf( "\n %s\n", desc );
  for( i = 0; i < m; i++ ) {
    for( j = 0; j < n; j++ ) printf( " %9.5f", a[i+j*lda] );
    printf( "\n" );
  }
}

float estimate_Wiener(fftwf_complex ** dataF, int nx, int ny, int z, int nphases, float rdistcutoff)
/*
To estimate signal to noise ratio, which is 1/wiener, based on the raw data's power
outside and inside of the resolution limit.
dataF is in right half Fourier space.
*/
{
  int ph, i, j, kx, ky, nxy, xdim;
  float kr;
  fftwf_complex * pIm;
  double sumSignalSquare=0, sumNoiseVariance=0;
  int countSignal=0, countNoise=0;

  xdim = (nx+2)/2;
  nxy = xdim*ny;

  /* for each of the phases in the raw data, estimate noise variance*/
  for (ph=0; ph<nphases; ph++) {
    pIm = dataF[ph] + z*nxy;
    for (i=0; i<ny; i++) {
      ky = i;
      if (ky>ny/2) ky -= ny;
      for (j=0; j<xdim; j++) {
        kx = j;
        kr = sqrt(kx*kx+ky*ky);
        if (kr<rdistcutoff) {
          countSignal ++;
          sumSignalSquare += mag2(pIm[j]);
        }
        else {
          countNoise ++;
          sumNoiseVariance += mag2(pIm[j]);
        }
      }
      /* advance pIm by one row */
      pIm += xdim;
    }
  }
  sumSignalSquare /= countSignal;
  sumNoiseVariance /= countNoise;

  return sumNoiseVariance / sumSignalSquare;
}

float estimate_Wiener_3D(fftwf_complex *dataF, int nx, int ny, int nz, float dkx, float dky, float dz, int wavelength, ReconParams* pParams)
/*
To estimate signal to noise ratio, which is 1/wiener, based on the raw data's power
outside the resolution limit and power at DC.
dataF is the real-to-complex FFT of band0
*/
{
  int i, j, k, nxy, xdim;
  float kr, kx, ky, dkz;
  fftwf_complex * pIm;
  double sumNoiseVariance = 0, sumMidFrequencies=0;
  int countNoise=0, countMidFrequencies=0;
  float lambdaem, alpha, rdistcutoff;
  int zdistcutoff;

  xdim = (nx+2)/2;
  nxy = xdim*ny;

  /* Calculate x-y and z resolution limit */
  
  dkz = 1./(nz*dz);
  lambdaem = (wavelength/pParams->nimm)/1000.0;  /* emission wavelength in the sample, in microns */
  alpha = asin(pParams->na/pParams->nimm);  /* aperture angle of objectives */
  rdistcutoff = (pParams->na*2/(wavelength/1000.0));
  zdistcutoff = (int) ceil(((1-cos(alpha))/lambdaem) / dkz);    /* OTF support axial limit in data pixels */

  for (k=zdistcutoff; k<nz-zdistcutoff; k++) {
    pIm = dataF + k*nxy;
    for (i=0; i<ny; i++) {
      ky = i;
      if (ky>ny/2) ky -= ny;
      ky *= dky;
      for (j=0; j<xdim; j++) {
        kx = j * dkx;
        kr = sqrt(kx*kx + ky*ky);
        if ( (k>zdistcutoff && k<nz-zdistcutoff) || kr>rdistcutoff) {
          countNoise ++;
          sumNoiseVariance += mag2(pIm[j]);
        }
        if ( (k<zdistcutoff || k>nz-zdistcutoff) && fabs(kr-rdistcutoff/2) < 0.03 ) { /* statistics of mid-frequency signal strength*/
          countMidFrequencies ++;
          sumMidFrequencies += mag2(pIm[j]);
        }
      }
      /* advance pIm by one row */
      pIm += xdim;
    }
  }

  sumNoiseVariance /= countNoise;
  printf("Normalized DC amplitude: %g\n", cabs(dataF[0])/(nx*ny*nz));

  printf("noise stdev: %g\n", sqrt(sumNoiseVariance));
  printf("Mid-frequency average: %g\n", sqrt(sumMidFrequencies/countMidFrequencies));

  return sqrt(sumNoiseVariance) / sqrt(sumMidFrequencies/countMidFrequencies) * 7.0e-3;
}

void calculateSectorBoundaries(int ndirs, vector* k0, float * boundaryAngles)
{
  int i;

  for (i=0; i<ndirs; i++) {

    float k0angle;

    k0angle = atan2(k0[i].y, k0[i].x);
    boundaryAngles[i] = k0angle - M_PI/ndirs/2;
    boundaryAngles[i+ndirs] = boundaryAngles[i] + M_PI;

    // atan2() returns between [-pi, pi]; so wrap the boundary angles inside this range as well
    if (boundaryAngles[i] > M_PI)
      boundaryAngles[i] -= 2*M_PI;
    else if (boundaryAngles[i] < -M_PI)
      boundaryAngles[i] += 2*M_PI;
    if (boundaryAngles[i+ndirs] > M_PI)
      boundaryAngles[i+ndirs] -= 2*M_PI;
    else if (boundaryAngles[i+ndirs] < -M_PI)
      boundaryAngles[i+ndirs] += 2*M_PI;
  }
  for (i=0; i<2*ndirs; i++)
    printf("%.3f   ", boundaryAngles[i]);
  printf("\n");
}


int angleIsInBetween(float angle, float angle0, float angle1)
// all angle values should be between [-pi, pi]
{
  float smaller, larger;
  int ret = 0;

  smaller = angle0 > angle1 ? angle1 : angle0;
  larger  = angle0 > angle1 ? angle0 : angle1;

  if (larger - smaller < M_PI/2) { // normal logic will be applied
    if (angle > smaller && angle <= larger)
      ret = 1;
  }
  else { // apply the opposite logic
    if (angle > 0 && angle > larger && angle-2*M_PI <= smaller ||
        angle < 0 && angle <= smaller && angle+2*M_PI > larger)
      ret = 1;
  }
  return ret;
}

void intersectLineCircle(float angle, float x0, float y0, float R, vector *intersection)
// find where a through-origin line, whose angle is given by angle, intersects a circle;
// (x0, y0) is the center of the circle of radius R
{
  float k = tan(angle);
  float numerator_part1 = x0 + k*y0;
  float numerator_part2 = sqrt(((x0+k*y0)*(x0+k*y0) - (1+k*k)*(x0*x0+y0*y0-R*R)));
  float denom = 1+k*k;

  float xsolution1 = (numerator_part1 + numerator_part2) / denom;
  float xsolution2 = (numerator_part1 - numerator_part2) / denom;

  intersection[0].x = xsolution1;
  intersection[0].y = xsolution1 * k;
  intersection[1].x = xsolution2;
  intersection[1].y = xsolution2 * k;
}

float getNextSectorStartingAngle(float *boundaryAngles, int dir, int ndirs, float deltaAngle)
{
  float roughTargetAngle;
  float aLargeRadian = 2*M_PI;
  int i, savedInd = dir;

  roughTargetAngle = boundaryAngles[dir] + deltaAngle;

  for (i=0; i<2*ndirs; i++) {
    float diffAngle = roughTargetAngle - boundaryAngles[i];
    if (diffAngle > M_PI)
      diffAngle -= 2*M_PI;
    else if (diffAngle < -M_PI)
      diffAngle += 2*M_PI;
    if (fabs(diffAngle) < aLargeRadian) {
      aLargeRadian = fabs(diffAngle);
      savedInd = i;
    }
  }
  return boundaryAngles[savedInd];
}


float getPreciseResLimit(float arg, float * boundaryAngles, vector *k0, int ndirs, int norders, float R)
// boundaryAngles has 2*ndirs elements, to take into account both k0 and -k0 scenarios
{
  int i;
  vector intersection[2];
  float dist[2];
  float deltaAngle = M_PI/ndirs;


  for (i=0; i<ndirs; i++) {
    float higher_angle, higher_angle_plus_PI;
    higher_angle = getNextSectorStartingAngle(boundaryAngles, i, ndirs, deltaAngle);
    higher_angle_plus_PI = higher_angle + M_PI;
    if (higher_angle_plus_PI > M_PI) 
      higher_angle_plus_PI -= 2*M_PI;
    else if (higher_angle_plus_PI < -M_PI)
      higher_angle_plus_PI += 2*M_PI;


    if (angleIsInBetween(arg, boundaryAngles[i], higher_angle) ||
        angleIsInBetween(arg, boundaryAngles[i+ndirs], higher_angle_plus_PI)) {
      intersectLineCircle(arg, k0[i].x * (norders-1), k0[i].y * (norders-1), R, intersection);
      break;
    }
  }
  for (i=0; i<2; i++)
    dist[i] = sqrt(intersection[i].x*intersection[i].x + intersection[i].y*intersection[i].y);

  return ( dist[0] > dist[1] ? dist[0] : dist[1] );
}
