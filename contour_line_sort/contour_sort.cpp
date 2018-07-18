/*
  This code rearranges the contour data exported from Paraview.
      -     finds the nearest-distance coordinate from a reference coordinate.

  Dong-Uk Kim, July 3rd 2018
*/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
  int nline = 1305;
  int ndat = nline;
  char datfname[100];
  char outfilename[200];
  sprintf(datfname, "Contour_as_exported.csv"); //input data
  sprintf(outfilename, "phase_field_1500_1.csv");  //output data
  double xinit = 0.0;
  double yinit = 400.0;
  double zinit = 0.0;
  ////////////////////
  FILE *fin;
  double *xin, *yin, *zin, *d1, *d2, *d3;
  double *xout, *yout;
  double *rNNd;
  double *sNNd;
  xin = new double [ndat];
  yin = new double [ndat];
  zin = new double [ndat];
  d1 = new double [ndat];
  d2 = new double [ndat];
  d3 = new double [ndat];
  xout = new double [ndat];
  yout = new double [ndat];
  rNNd = new double [ndat];
  sNNd = new double [ndat];
  ////////////////////
  char linebuff[100];
  ////////////////////
  void print_3_darrays(double *a1, double *a2, double *a3, int nsize);
  void pick_mindist_coord(double *xout, double *yout, double xref, double yref, double *xsample, double *ysample, int iref, int istart, int iend);
  void calc_NNdist_arr(double *output, double *xarr, double *yarr, int nsize);
  void fileout(char *fname, double *x, double *y, int ilb, int iub);
  ////////////////////

  //data loading and indexing (in order of presence in the file)
  fin = fopen(datfname,"r");
  fscanf(fin,"%s", linebuff);
  for (int i=1; i<ndat; i++) {  //The 0th index is reserved for the initial reference coordinate that user will enter
    fscanf(fin, "%lf,%lf,%lf,%lf,%lf,%lf",&xin[i],&yin[i],&zin[i],&d1[i],&d2[i],&d3[i]);
  }
  fclose(fin);

  //Saving the initial reference coordinate
  xin[0] = xinit;
  yin[0] = yinit;
  zin[0] = zinit;

  //Printing the input data on the screen
  print_3_darrays(xin,yin,zin,ndat);

  //initialization of xinter and yinter for the forward swipping
  xout[0]=xin[0];
  yout[0]=yin[0];
  
  //Find the nearest coordinate from i-1 coordinate, forward swipping
  for (int i = 1; i<ndat; i++) {
    double xmin, ymin;
    int iref = i-1;
    int iscan_lb = i;
    int iscan_ub = ndat-1;
    pick_mindist_coord(&xmin, &ymin, xin[iref], yin[iref], xin, yin, iref, iscan_lb, iscan_ub);
    xout[i] = xmin;
    yout[i] = ymin;
  }
  calc_NNdist_arr(sNNd, xout, yout, ndat);
  print_3_darrays(xout,yout,sNNd,ndat);

  fileout(outfilename, xout,yout,1,ndat-1);
  

  return 0;
}

void print_3_darrays(double *a1, double *a2, double *a3, int nsize)
{
  for (int i=0; i<nsize; i++){
    printf("%6d %15lf %15lf %15lf\n",i, a1[i], a2[i], a3[i]);
  }
}

void calc_NNdist_arr(double *output, double *xarr, double *yarr, int nsize)
{
  double dist(double x1, double y1, double x2, double y2);
  output[0] = 0.;
  for (int i=0+1; i<nsize; i++){
    //output[i] = sqrt(pow(xref-xarr[i],2)+pow(yref-yarr[i],2));
    output[i] = dist(xarr[i-1],yarr[i-1],xarr[i],yarr[i]);
  }
}

double dist(double x1, double y1, double x2, double y2)
{
  return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}

void pick_mindist_coord(double *xout, double *yout, double xref, double yref, double *xsample, double *ysample, int iref, int istart, int iend)
{
  //////////////////////////
  int nitem = iend-istart+1;
  //////////////////////////
  double *d;
  d = new double [nitem];
  //////////////////////////
  double dmin=3.0857e16; //1 Parsec in SI unit; the largest length constant I know.
  double xmin, ymin;
  int imin;
  //////////////////////////
  double dist(double x1, double y1, double x2, double y2);
  void switching_data(double *x, double *y, int i1, int i2);
  //////////////////////////
  //make a distance list
  for (int i=istart, id=0; i<=iend; i++, id++){
    d[id] = dist(xref,yref, xsample[i], ysample[i]);
    if (d[id] <= dmin) {
      dmin = d[id];
      xmin = xsample[i];
      ymin = ysample[i];
      imin = i;
    }
  }
  //switching iref data with imin data
  switching_data(xsample,ysample, istart, imin);
  *xout = xmin;
  *yout = ymin;
}

void switching_data(double *x, double *y, int i1, int i2)
{
  /////////////////////////
  int itmp;
  double xtmp;
  double ytmp;
  /////////////////////////
  xtmp = x[i1];
  ytmp = y[i1];
  x[i1] = x[i2];
  y[i1] = y[i2];
  x[i2] = xtmp;
  y[i2] = ytmp;
}

void fileout(char fname[200], double *x, double *y, int ilb, int iub)
{
  FILE *fout;
  
  fout = fopen(fname, "w");
  fprintf(fout,"x,y\n");
  for (int i=ilb; i<=iub; i++) {
    fprintf(fout,"%lf,%lf\n",x[i],y[i]);
  }
  fclose(fout);
}
