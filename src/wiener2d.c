#include <stdio.h>
#include <stdlib.h>  /* free, malloc, strtol */
#include <string.h>
#include <unistd.h>  /* getopt */
#include <errno.h>   /* errno */
#include <IMInclude.h>
#include <math.h>
#include <sfftw.h>
#include <srfftw.h>

float wiener=1e-4;


int scan_float(const char *text, float *p_value)
{
  float scanned;

  errno = 0;  /* defined in errno.h */
  scanned = strtod(text, NULL);
  if (errno ==0)
	*p_value = scanned;
  else return 1;

  return 0;
}

int scan_int(const char *text, int *p_value)
{
  int scanned;

  errno = 0;  /* defined in errno.h */
  scanned = strtol(text, NULL, 10);
  if (errno ==0)
	*p_value = scanned;
  else return 1;

  return 0;
}

fftw_complex conjugate(fftw_complex a)
{
  a.im *= -1;
  return a;
}

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
  int iin, jin, indin, kx, ky, kycent, kxcent;
  fftw_complex exp_iphi;
  float phi1, phi, dphiy, dphix;

  kycent = ny/2;
  kxcent = nx/2;

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
  int i, j, maxi=-1, maxj=-1, ind;
  int iminus, iplus, jminus, jplus;
  float maxval, reval, valminus, valplus;
  double sum;

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
  for (i=0; i<*yc-20; i++)
	for (j=0; j<nx; j++)
	  sum += rawimage[i*nx + j];
  *background = sum / ((*yc-20)*nx);
}

void fft_rescale(float *p, int nx, int ny, int nz, int nx_plus_2)
{
  float scale=1.0/(nx*ny*nz);
  int actual_nxy, actual_nx, k, i, j, rowstart, secstart;

  if (nx_plus_2)
    actual_nx = nx+2;
  else
    actual_nx = nx;
  actual_nxy = actual_nx*ny;

  for (k=0; k<nz; k++) {
    secstart = k * actual_nxy;
    for (i=0; i<ny; i++) {
      rowstart = secstart + i*actual_nx;
      for (j=0; j<nx; j++)
	p[rowstart+j] *= scale;
    }
  }
}


void fixorigin2D(fftw_complex * otfkx, int kx1, int kx2)
{
  float meani, slope, avg, lineval;
  int numvals, i;
  double *sum, totsum=0, ysum=0, sqsum=0;

  meani = 0.5*(kx1+kx2);
  numvals = kx2-kx1+1;
  sum = (double *) malloc((kx2+1)*sizeof(double));

  for (i=0; i<=kx2; i++) {  /* '<=' is correct. For a long time, '<' was used, which is a bug. */
	sum[i] = cmag(otfkx[i]);
    if (i>=kx1) {
      totsum += sum[i];
      ysum += sum[i] * (i-meani);
      sqsum += (i-meani) * (i-meani);
    }
  }
  slope = ysum / sqsum;
  avg = totsum / numvals;

  for (i=0; i<kx1; i++) {
    lineval = avg + (i-meani)*slope;
    otfkx[i].re -= sum[i] - lineval;
  }

  free(sum);
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
		avg_output[indout].re += band[indin].re;
		avg_output[indout].im += band[indin].im;
		count[indout] ++;
	  }
	}
  }
  
  for (indout=0; indout<nx/2+1; indout++) {
    if (count[indout]>0) {
      avg_output[indout].re /= count[indout];
      avg_output[indout].im /= count[indout];
    }
  }

  free(count);
}

fftw_complex otfinterpolate(fftw_complex * otf, float rindex, int nrotf)
{
  int irindex, index00, index01;
  fftw_complex otfval;
  float ar;

  if (rindex < 0 ) {
	printf("rindex is less than 0. rindex=%f\n", rindex);
	exit(-1);
  }

  irindex = floor(rindex);

  ar = rindex - irindex;

  index00 = irindex;
  index01 = irindex+1;

  otfval.re = otf[index00].re*(1-ar) + otf[index01].re*ar;
  otfval.im = otf[index00].im*(1-ar) + otf[index01].im*ar;
    
  return otfval;
}


void wienerfilter(fftw_complex * g, int nx, int ny, float dkr, fftw_complex * otf, int nr_otf,
				  float dkr_otf, float rcutoff, fftw_complex * output)
{
  /* 'g' is the raw data's FFT */
  int i, j, ind, irel;
  float amp2, rho, kr, krscale;
  fftw_complex A_star_g, otf_val;
  float w;

  memset(output, 0, (nx/2+1)*ny*sizeof(fftw_complex));
  
  w = wiener*wiener;
  krscale = dkr/dkr_otf;

  for (i=0; i<ny; i++) {
    irel = i;
    if (irel > ny/2)
      irel -= ny;
    for (j=0; j<nx/2+1; j++) {
      ind = i*(nx/2+1) + j;
	  kr = sqrt(irel*irel+j*j);
	  if (kr <=rcutoff) {
		otf_val = otfinterpolate(otf, kr*krscale, nr_otf);

		amp2 = otf_val.re*otf_val.re + otf_val.im*otf_val.im;
		A_star_g = cmult(conjugate(otf_val), g[ind]);

		/* apodization */
		rho = sqrt(irel*irel+j*j)/rcutoff;

		output[ind].re = A_star_g.re / (amp2+w) * (1-rho);
		output[ind].im = A_star_g.im / (amp2+w) * (1-rho);
	  }
	}
  }
}

void apodize(int napodize,int nx,int ny,float *image)
{
  float diff,fact;
  int k,l;

  for(k=0;k<nx;k++)	{
	diff = (image[(ny-1)*nx +k] - image[ /* 0*nx+ */ k ]) / 2;
	for(l=0; l<napodize; l++) {
	  fact = 1 - sin((((float)l+0.5)/napodize)*M_PI*0.5);
	  image[l*nx + k] += diff*fact;
	  image[(ny-1-l)*nx + k] -= diff*fact;
	}
  }
  for(l=0;l<ny;l++)	{
	diff = (image[l*nx + nx-1] - image[l*nx /* +0 */]) / 2;
	for(k=0; k<napodize; k++) {
	  fact = 1 - sin((((float)k+0.5)/napodize)*M_PI*0.5);
	  image[l*nx + k] += diff*fact;
	  image[l*nx + ny-1-k] -= diff*fact; 
	}
  }
}

float average(float *buf, int size)
{
  float sum=0, *p=buf;
  int i;

  for (i=0; i<size; i++)
    sum += *p++;
  return sum / size;
}


int main(int argc, char *argv[])
{
  int istream_no, ostream_no, psfstream_no, avgstream_no;
  int i, j, t, nx_psf, nx_otf, ny_psf, nx, ny, nz, nt, ixyz[3], mxyz[3], pixeltype, napodize = 10;
  float xcofm, ycofm, background_psf, background=200, min, max, mean, dr, dr_psf, dkr, dkr_otf;
  float *psf_array, *raw_image, *proc_image, *buffer=NULL;
  fftw_complex * otf2d, *raw_fft, *proc_fft, *otf_avg;
  rfftwnd_plan rfftplan, rfftplan_inv;
  IW_MRC_HEADER header, psf_header;
  int nxy;
  float NA=1.2, rdistcutoff;
  int nphases = 1;  /* allow input to be 2D SIM datasets */
  float ref=0; /* used in fading correction for time-series */

  extern char *optarg;
  char letter, otfoutput[120];

  otfoutput[0]='\0';  /* empty */

  IMAlPrt(0);       /* suppress printout of file header info */

  while ((letter=getopt(argc,argv,"w:b:hO:N:p:"))!= -1) {
	if (letter == 'w') {
	  if (scan_float(optarg, &wiener) != 0 || wiener <= 0) {
		fprintf(stderr, "wiener constant has to be positive\n");
		printf("%f\n", wiener);
		return 1;
	  }
	}
	else if (letter == 'h') {
	  fprintf(stderr, "Usage: wiener2d input_file output_file psf_file [-w wiener] [-b background] [-N NA] [-O OTF2save] [-h]\n");
	  return 0;
	}
	else if (letter == 'b') {
	  if (scan_float(optarg, &background) != 0) {
		fprintf(stderr, "invalid argument supplied to -b option\n");
		return 1;
	  }
	}
	else if (letter == 'N') {
	  if (scan_float(optarg, &NA) != 0) {
		fprintf(stderr, "invalid argument supplied to -N option\n");
		return 1;
	  }
	}
	else if (letter == 'O') {
	  if (strcpy(otfoutput, optarg) != otfoutput) {
		fprintf(stderr, "invalid argument supplied to -O option\n");
		return 1;
	  }
    }
	else if (letter == 'p') {
	  if (scan_int(optarg, &nphases) != 0) {
		fprintf(stderr, "invalid argument supplied to -p option\n");
		return 1;
	  }
	}
  }

  if (argc < 4) {
	fprintf(stderr, "Usage: wiener2d input_file output_file psf_file [-w wiener] [-b background] [-N NA] [-p nphases] [-O OTF2save] [-h]\n");
	return 0;
  }

  istream_no = 1;
  ostream_no = 2;
  psfstream_no = 3;
  avgstream_no = 4;

  if (IMOpen(istream_no, argv[optind], "ro")) {
	printf("Reading input file %s\n", argv[optind]);
    fprintf(stderr, "File %s does not exist.\n", argv[optind]);
    return -1;
  }

  if (IMOpen(ostream_no, argv[optind+1], "new")) {
	printf("Creating output file %s\n", argv[optind+1]);
    fprintf(stderr, "Couldn't write to file %s.\n", argv[optind+1]);
    return -1;
  }

  if (IMOpen(avgstream_no, "averaged.mrc", "new")) {
    fprintf(stderr, "Couldn't write to average file.\n");
    return -1;
  }

  if (IMOpen(psfstream_no, argv[optind+2], "ro")) {
	printf("Reading psf file %s\n", argv[optind+2]);
    fprintf(stderr, "File %s does not exist.\n", argv[optind+2]);
    return -1;
  }

  printf("wiener = %.3e\n", wiener);
  printf("background = %.1f\n", background);

  IMRdHdr(istream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(istream_no, &header);

  nx = header.nx;
  ny = header.ny;
  nz = header.nz;
  nt = header.num_times;
  dr = header.xlen;
  nz /= nphases*nt;  /* nz is more like "ndirs" in SIM */
  printf("nz=%d\n", nz);

  nxy = nx * ny;

  IMRdHdr(psfstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(psfstream_no, &psf_header);

  if (psf_header.mode != IW_COMPLEX) {
    nx_psf = psf_header.nx;
    ny_psf = psf_header.ny;
    dr_psf = psf_header.xlen;

    nx_otf = nx_psf/2+1;

    /* Now load PSF data */
    psf_array = (float *) calloc(nx_psf*ny_psf, sizeof(float));
    IMRdSec(psfstream_no, psf_array);

    determine_center_and_background(psf_array, nx_psf, ny_psf, &xcofm, &ycofm, &background_psf);

    printf("PSF's center of mass is (%.3f, %.3f)\n Background is %.3f\n", xcofm, ycofm, background_psf);

    /* subtract background from PSF */
    for (i=0; i<nx_psf*ny_psf; i++)
      psf_array[i] -= background_psf;

    otf2d = calloc(nx_otf*ny_psf, sizeof(fftw_complex));

    rfftplan = rfftw2d_create_plan(ny_psf, nx_psf, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
    rfftwnd_one_real_to_complex(rfftplan, psf_array, otf2d);
    fftwnd_destroy_plan(rfftplan);
  
    shift_center(otf2d, nx_psf, ny_psf, xcofm, ycofm);
  
    otf_avg = calloc(nx_otf, sizeof(fftw_complex));
    radialft(otf2d, nx_psf, ny_psf,otf_avg);

    /* normalize OTF */
    dkr_otf = 1.0/(dr_psf*nx_psf);

    /* clean up outside the observable region */
    rdistcutoff = 2*NA/(psf_header.iwav1*0.001) / dkr_otf;
    printf("dr_psf=%.3f, nx_psf=%d, wavelength= %d, rdistcutoff = %.3f\n", dr_psf, nx_psf, psf_header.iwav1, rdistcutoff);

    if (rdistcutoff > nx_psf/2) {
      char s[100];
      fprintf(stderr, "Data undersampled. What cutoff value should be used? ");
      fgets(s, 100, stdin);
      sscanf(s, "%f", &rdistcutoff);
      printf("Using user-provided rdistcutoff: %.3f\n", rdistcutoff);
    }

    for (i=(int)rint(rdistcutoff); i<nx_otf; i++) {
      otf_avg[i].re = 0;
      otf_avg[i].im = 0;
    }

    fixorigin2D(otf_avg, 1, 12);

    for (i=1; i<nx_otf; i++) {
      otf_avg[i].re /= otf_avg[0].re;
      otf_avg[i].im /= otf_avg[0].re;
    }
    otf_avg[0].re = 1; otf_avg[0].im = 0;

    if (otfoutput[0] != '\0') {
      FILE *f=fopen(otfoutput, "wb");
      fwrite(otf_avg, sizeof(fftw_complex), nx_otf, f);
      fclose(f);
      f=fopen(strcat(otfoutput, ".beforeRadialft"), "wb");
      fwrite(otf2d, sizeof(fftw_complex), ny_psf*nx_otf, f);
      fclose(f);
    }

    free(otf2d);
    free(psf_array);
  }
  else { /* (psf_header.mode == IW_COMPLEX) */
    nx_otf = psf_header.nx;
    dkr_otf = psf_header.xlen;
    otf_avg = calloc(nx_otf, sizeof(fftw_complex));
    IMRdSec(psfstream_no, otf_avg);
  }
  
  IMClose(psfstream_no);

  raw_image = malloc(nxy*sizeof(float));
  if (nphases > 1)
    buffer = calloc(nxy, sizeof(float));
  raw_fft = malloc((nx+2)/2*ny*sizeof(fftw_complex));

  rfftplan = rfftw2d_create_plan(ny, nx, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  rfftplan_inv = rfftw2d_create_plan(ny, nx, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);

  proc_image = calloc(nxy, sizeof(float));
  proc_fft = calloc((nx+2)/2*ny, sizeof(fftw_complex));
  dkr = 1.0/(dr*nx);
  rdistcutoff = 2*NA/(header.iwav1*0.001) / dkr;


  IMTrHdr(ostream_no, istream_no);
  IMAlMode(ostream_no, IW_FLOAT);
  ixyz[0] = nx;
  ixyz[1] = ny;
  /* if (nphases>1) */
  /*   ixyz[2] = (nz+1)*nt; */
  /* else */
  ixyz[2] = nz*nt;
  IMAlSiz(ostream_no, ixyz, mxyz);

  IMTrHdr(avgstream_no, istream_no);
  IMAlMode(avgstream_no, IW_FLOAT);
  ixyz[2] = nt;
  IMAlSiz(avgstream_no, ixyz, mxyz);

  for (t=0; t<nt; t++)
    for (i=0; i<nz; i++) {  /* when nphases > 1, nz means "ndirs" */
      if (nphases ==1) {
        IMRdSec(istream_no, raw_image);
        /* subtract background */
        for (j=0; j<nxy; j++)
          raw_image[j] -= background;
      }
      else {
        int p;
        memset(raw_image, 0, nxy*sizeof(float));
        /* Input is SIM dataset */
        for (p=0; p<nphases; p++) {
          IMRdSec(istream_no, buffer);
          for (j=0; j<nxy; j++)
            raw_image[j] += (buffer[j] - background)/nphases;
        }
      }

      if (t==0 && i==0)
        ref = average(raw_image, nxy);
      else {
        float ratio;
        ratio = ref / average(raw_image, nxy);
        for (j=0; j<nxy; j++)
          raw_image[j] *= ratio;
      }

      if (nphases > 1 && i==0)
        IMWrSec(avgstream_no, raw_image);  // save the averaged image in a separate output file

        
      apodize(napodize, nx, ny, raw_image);

      rfftwnd_one_real_to_complex(rfftplan, raw_image, raw_fft);
      wienerfilter(raw_fft, nx, ny, dkr, otf_avg, nx_otf, dkr_otf, rdistcutoff, proc_fft);
      rfftwnd_one_complex_to_real(rfftplan_inv, proc_fft, proc_image);
      fft_rescale(proc_image, nx, ny, 1, 0);

      printf("section %d done.\n", i);
      IMWrSec(ostream_no, proc_image);
    }

  fftwnd_destroy_plan(rfftplan);
  fftwnd_destroy_plan(rfftplan_inv);
  
  IMClose(istream_no);

  IMWrHdr(ostream_no, "Processed", 1, min, 1000, mean);
  IMClose(ostream_no);
  
  IMWrHdr(avgstream_no, "averaged", 1, min, 1000, mean);
  IMClose(avgstream_no);
  
  return 0; 
}
