#include <stdio.h>
#include <stdlib.h>  /* free, malloc, strtol */
#include <string.h>
#include <unistd.h>  /* getopt */
#include <errno.h>   /* errno */
#include <IMInclude.h>
#include <math.h>
#include <sfftw.h>
#include <srfftw.h>
#include <tiffio.h>


fftw_complex cmult(fftw_complex a, fftw_complex b)
{
  fftw_complex c;
  c.re = a.re*b.re - a.im*b.im;
  c.im = a.re*b.im + a.im*b.re;
  return(c); 
}

float cmag(fftw_complex a)
{
  return sqrt(a.re*a.re + a.im*a.im);
}


void shift_center(fftw_complex *bands, int nx, int ny, float xc, float yc)
{
  int iin, jin, indin, nxy, kx, ky, kycent, kxcent;
  fftw_complex exp_iphi;
  float phi1, phi, dphiy, dphix;

  kycent = ny/2;
  kxcent = nx/2;
  nxy = (nx/2+1)*ny;

  dphiy = 2*M_PI*yc/ny;
  dphix = 2*M_PI*xc/nx;

  for (iin=0; iin<ny; iin++) {
    ky = iin;
    if (iin>kycent) ky -= ny;
    phi1 = dphiy*ky;   /* first part of phi */
    for (jin=0; jin<kxcent+1; jin++) {
      kx = jin;
      indin = iin*(nx/2+1)+jin;
      phi = phi1+dphix*kx;  /* second part of phi */

      exp_iphi.re = cos(phi);
      exp_iphi.im = sin(phi);
      bands[indin] = cmult(bands[indin], exp_iphi);
    }
  }
}

float fitparabola( float a1, float a2, float a3 )
{
  float slope,curve,peak;

  slope = 0.5* (a3-a1);         /* the slope at (x=0). */
  curve = (a3+a1) - 2*a2;       /* (a3-a2)-(a2-a1). The change in slope per unit of x. */
  if( curve == 0 ) 
    {
      printf("no peak: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f\n",a1,a2,a3,slope,curve);
      return( 0.0 );
    } 
  peak = -slope/curve;          /* the x value where slope = 0  */
  if( peak>1.5 || peak<-1.5 )
    {
      printf("bad peak position: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f, peak=%f\n",a1,a2,a3,slope,curve,peak);
      return( 0.0 );
    } 
  return( peak );
}

/*  locate peak pixel to subpixel accuracy by fitting parabolas  */
void determine_center_and_background(float *rawimage, int nx, int ny, float *xc, float *yc, float *background)
{
  int i, j, maxi=0, maxj=0, ind, nxy2;
  int iminus, iplus, jminus, jplus;
  float maxval, reval, valminus, valplus, minval;
  double sum;

  /*   printf("In determine_center_and_background()\n"); */
  nxy2 = nx*ny;

  /* Search for the peak pixel */

  maxval=0.0;
  for(i=0;i<ny;i++)
    for(j=0;j<nx;j++) {
      ind=i*nx+j;
      reval=rawimage[ind];
      if( reval > maxval ) {
        maxval = reval;
        maxi=i; maxj=j;
      }
    }

  /*   printf ("maxi=%d, maxj=%d, maxval=%f", maxi, maxj, maxval); */

  iminus = maxi-1; iplus = maxi+1;
  if( iminus<0 ) iminus+=ny;
  if( iplus>=ny ) iplus-=ny;
  jminus = maxj-1; jplus = maxj+1;
  if( jminus<0 ) jminus+=nx;
  if( jplus>=nx ) jplus-=nx;

  valminus = rawimage[iminus*nx+maxj];
  valplus  = rawimage[iplus *nx+maxj];
  *yc = maxi + fitparabola(valminus, maxval, valplus);

  valminus = rawimage[maxi*nx+jminus];
  valplus  = rawimage[maxi*nx+jplus];
  *xc = maxj + fitparabola(valminus, maxval, valplus);
  
  sum = 0;
  minval = maxval;
  for (i=0; i<maxi-20; i++)
    for (j=0; j<nx; j++)
      sum += rawimage[i*nx + j];
  /*      if (rawimage[i*nx+j]<minval) */
  /*        minval = rawimage[i*nx+j]; */
  *background = sum / ((maxi-20)*nx);
  /*   *background = minval; */
}

void fixorigin2D(fftw_complex * otfkx, int kx1, int kx2)
{
  float meani, slope, avg, lineval;
  int numvals, i;
  double  val, totsum=0, ysum=0, sqsum=0;

  meani = 0.5*(kx1+kx2);
  numvals = kx2-kx1+1;

  for (i=kx1; i<=kx2; i++) {
    val = otfkx[i].re;
    totsum += val;
    ysum += val * (i-meani);
    sqsum += (i-meani) * (i-meani);
  }
  slope = ysum / sqsum;
  avg = totsum / numvals;

  for (i=0; i<kx1; i++) {
    lineval = avg + (i-meani)*slope;
    otfkx[i].re = lineval;
    otfkx[i].im = 0;
  }

}

void radialft(fftw_complex *band, int nx, int ny, fftw_complex *avg_output)
{
  int iin, jin, indin, indout, kx, ky, kycent, kxcent;
  int *count;
  float rdist;

  /*   printf("In radialft()\n"); */
  kycent = ny/2;
  kxcent = nx/2;

  if (!(count = (int *) calloc(nx/2+1, sizeof(int)))) {
    fprintf(stderr, "Cannot allocate memory for count in radialft()\n");
    exit(-1);
  }

  for (iin=0; iin<ny; iin++) {
    ky = iin;
    if (iin>kycent) ky -= ny;
    for (jin=0; jin<kxcent+1; jin++) {
      kx = jin;
      rdist = sqrt(kx*kx+ky*ky);
      if (rint(rdist) < nx/2+1) {
        indin = iin*(nx/2+1)+jin;
        indout = rint(rdist);

        if (jin == 0) {
          avg_output[indout].re += band[indin].re;
          count[indout] ++;
        } else {     /* enter the effect of the diametrically opposite pixel, which contains the complex conjugate value */
          avg_output[indout].re += 2 * band[indin].re;
          count[indout] += 2;
        }
        /* avg_output[indout].im += band[indin].im;*/
        avg_output[indout].im = 0;  /* averaged im part must be zero by symmetry: im parts of opposite pixels have
                                       equal magnitudes and opposite signs*/

      }
    }
  }
  
  for (indout=0; indout<nx/2+1; indout++) {
    if (count[indout]>0) {
      avg_output[indout].re /= count[indout];
      /* avg_output[indout].im /= count[indout]; */   /* im part is already zero; does not need to be normalized */
    }
  }

  free(count);
}

int main(int argc, char *argv[])
{
  int psfstream_no, ostream_no, ixyz[3], mxyz[3], pixeltype, nx, ny, nz, i;
  IW_MRC_HEADER header, otf_header;
  float *psf_array, dr, xcofm, ycofm, background_psf, min, max, mean, ref;
  fftw_complex * otf2d, *otf2d_avg=0;
  rfftwnd_plan rfftplan;
  int saving2D = 0; /* by default, save 1D averaged OTF */
  char letter;
  int krL, krH; /* the kr range between [krL, krH] that's used in fixorigin()*/
  float NA, dkr, krcutoff, lambda;

  TIFF *input_tiff; //, *output_tiff, *otf_tiff;
  int bTIFF = 0;

  IMAlPrt(0);

  NA=1.2;
  krL = 3;
  krH = 10;
  while ((letter=getopt(argc,argv,"ht2L:H:N:d:w:"))!= -1) {
    if (letter == '2') {
      saving2D = 1;
    }
    else if (letter == 'L') {
      /* lower kr value of the kr range used in fixorigin() */
      krL = strtol(optarg, NULL, 10);
    }
    else if (letter == 'H') {
      /* higher kr value of the kr range used in fixorigin() */
      krH = strtol(optarg, NULL, 10);
    }
    else if (letter == 'N') {
      /* detection NA */
      NA = strtof(optarg, NULL);
    }
    else if (letter == 't') {
      /* Input is TIFF */
      bTIFF = 1;
    }
    else if (letter == 'd') {
      /* pixel size in um (only valid for TIFF input) */
      dr = strtof(optarg, NULL);
    }
    else if (letter == 'w') {
      /* wavelength in um (only valid for TIFF input) */
      lambda = strtof(optarg, NULL);
    }
    else if (letter == 'h') {
      fprintf(stderr, "Usage: otf2d [-2] [-L smaller_kr] [-H greater_kr] [-N detection_NA] [-t] [-h] [-d pixel_size_in_um] [-w wavelength_in_um] PSF_input OTF_output\n");
      return 0;
    }

  }
  if (argc < 3) {
    fprintf(stderr, "Usage: otf2d  [-2] [-L smaller_kr] [-H greater_kr] [-N detection_NA] [-t] [-h] [-d pixel_size_in_um] [-w wavelength_in_um] PSF_input OTF_output\n");
    return 0;
  }

  if (!bTIFF) {
    psfstream_no = 1;

    if (IMOpen(psfstream_no, argv[optind], "ro")) {
      fprintf(stderr, "File %s does not exist.\n", argv[optind]);
      return -1;
    }


    IMRdHdr(psfstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
    IMGetHdr(psfstream_no, &header);

    nx = header.nx;
    ny = header.ny;
    nz = header.nz;
    dr = header.xlen;
    dkr = 1/(nx*dr);
    lambda = header.iwav1/1000.0;
  }
  else { /* handling TIFF input */
    if (!(input_tiff = TIFFOpen(argv[optind], "r"))) {
      fprintf(stderr, "Can't open input TIFF file %s.\n", argv[optind+1]);
      return -1;
    }
    TIFFGetField(input_tiff,TIFFTAG_IMAGEWIDTH,&nx);
    TIFFGetField(input_tiff,TIFFTAG_IMAGELENGTH,&ny);
    nz=0;
    do ++nz; while (TIFFReadDirectory(input_tiff));
    dkr = 1/(nx*dr);
  }

  if (nz>1)
    fprintf(stderr, "Warning: more than one 2D images in input.\n");
  
  psf_array = (float *) calloc(nx*ny, sizeof(float));

  if (!bTIFF) {
    IMRdSec(psfstream_no, psf_array);
    IMClose(psfstream_no);
  }
  else {
    load_tiff(input_tiff, 0, 0, psf_array);
    TIFFClose(input_tiff);
  }

  determine_center_and_background(psf_array, nx, ny, &xcofm, &ycofm, &background_psf);

  printf("PSF's center of mass is (%.3f, %.3f)\n Background is %.3f\n", xcofm, ycofm, background_psf);

  /* subtract background from PSF */
  for (i=0; i<nx*ny; i++)
    psf_array[i] -= background_psf;

  otf2d = calloc((nx+2)/2*ny, sizeof(fftw_complex));

  rfftplan = rfftw2d_create_plan(ny, nx, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  rfftwnd_one_real_to_complex(rfftplan, psf_array, otf2d);
  fftwnd_destroy_plan(rfftplan);
  shift_center(otf2d, nx, ny, xcofm, ycofm);


  dkr = 1/(nx*dr);
  krcutoff = 2*NA/lambda/dkr + .5;
  printf("NA=%f, krcutoff=%f\n", NA, krcutoff);

  if (saving2D) {
    ref = otf2d[0].re;
    for (i=0; i<(nx/2+1)*ny; i++) {
      otf2d[i].re /= ref;
      otf2d[i].im /= ref;
    } 
  } else {
    otf2d_avg = calloc(nx/2+1, sizeof(fftw_complex));
    radialft(otf2d, nx, ny, otf2d_avg);

    fixorigin2D(otf2d_avg, krL, krH);

    ref = otf2d_avg[0].re;
    for (i=0; i<krcutoff /*(nx/2+1)*/; i++) {
      otf2d_avg[i].re /= ref;
      otf2d_avg[i].im /= ref;
    }
    for (i=krcutoff; i<nx/2+1; i++) {
      otf2d_avg[i].re = 0;
      otf2d_avg[i].im = 0;
    }
    
  }
  
  ostream_no = 2;
  if (IMOpen(ostream_no, argv[optind+1], "new")) {
    fprintf(stderr, "Couldn't write to file %s.\n", argv[optind+1]);
    return -1;
  }
      
  if (!bTIFF) {
    otf_header = header;
    otf_header.nx = nx/2+1;
    if (saving2D) 
      otf_header.ny = ny;
    else
      otf_header.ny = 1;
    otf_header.xlen = 1/(nx*dr);
    otf_header.ylen = 1/(ny*dr);

    IMPutHdr(ostream_no, &otf_header);
    IMAlMode(ostream_no, IW_COMPLEX);
  }
  else {
    int ixyz[3];
    int mxyz[3] = {1, 1, 1};
    float pixelsize[3];
    float wavelengths[1];

    ixyz[0] = nx/2+1;
    if (saving2D)
      ixyz[1] = ny;
    else
      ixyz[1] = 1;
    ixyz[2] = 1;
    pixelsize[0] =  1/(nx*dr);
    pixelsize[1] = pixelsize[0];
    pixelsize[2] = 1;
    
    IMCrHdr(ostream_no, ixyz, mxyz, IW_COMPLEX, "", 1);
    IMAlDel(ostream_no, pixelsize);
    wavelengths[0] = (int) (lambda*1000);
    IMAlWav(ostream_no, 1, wavelengths);
  }

  if (saving2D)
    IMWrSec(ostream_no, otf2d);
  else
    IMWrSec(ostream_no, otf2d_avg);

  IMWrHdr(ostream_no, "Convereted from TIFF PSF file", 1, 0, 1, 0);
  IMClose(ostream_no);
  return 0;
}
