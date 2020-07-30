/***********************************************************/
/* All rights reserved by Amir Hossein Alikhah Mishamandani*/
/******* Created for ADCS system of SSS-P1 Satellite *******/
/*********************** 2018 - 2019 ***********************/
/************ World Magnetic Model 2015 - 2020 *************/
/*********** Original Program By: Dr. John Quinn ***********/
/********** Corrected by: Stefan Maus & Manoj Nair *********/
/***********************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "wmm.h"

#define NaN log(-1.0)

static int my_isnan(double d)
{
  return (d != d);              /* IEEE: only NaN is not equal to itself */
}

static void geomag(const char* cof_data, int *maxdeg);
static void geomg1(const char* cof_data, double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv);

void WorldMagneticModel(const char* cof_data, geomag_vector *inout, double dlat,double dlon,double altm,double time){
	memset(inout, 0, sizeof(geomag_vector));
	inout->X = 0;
	inout->Y = 0;
	inout->Z = 0;
	inout->F = 0;
	inout->H = 0;
	inout->decld = 0;
	inout->declm = 0;
	inout->incld = 0;
	inout->inclm = 0;	
	int   warn_H, warn_H_strong, warn_P;
	static int maxdeg;
	static double ati, adec, adip;
	static double alt, dec, dip, ti, gv;
	static double time1, dec1, dip1, ti1;
	static double time2, dec2, dip2, ti2;
	char decd[7], dipd[7],d_str[81],modl[20];
	char goodbye[81];
	double x1,x2,y1,y2,z1,z2,h1,h2;
	double ax,ay,az,ah;
	double rTd=0.017453292;
	double epochlowlim,epochuplim;
	double epochrange = 5.0;
	double warn_H_val, warn_H_strong_val;
	double dmin, imin, ddeg, ideg;

	if (sscanf(cof_data, "%lf%s",&epochlowlim,modl) < 2){
		fprintf(stderr, "Invalid header in model data.\n");
		exit(1);
	}
	
	maxdeg = 12;
    warn_H = 0;
    warn_H_val = 99999.0;
    warn_H_strong = 0;
    warn_H_strong_val = 99999.0;
    warn_P = 0;
    geomag(cof_data, &maxdeg);
	
    geomg1(cof_data, alt,dlat,dlon,time,&dec,&dip,&ti,&gv);
    time1 = time;
    dec1 = dec;
    dip1 = dip;
    ti1 = ti;
    time = time1 + 1.0;
      
    geomg1(cof_data, alt,dlat,dlon,time,&dec,&dip,&ti,&gv);
    time2 = time;
    dec2 = dec;
    dip2 = dip;
    ti2 = ti;
      
	/*COMPUTE X, Y, Z, AND H COMPONENTS OF THE MAGNETIC FIELD*/ 
    x1=ti1*(cos((dec1*rTd))*cos((dip1*rTd)));
    x2=ti2*(cos((dec2*rTd))*cos((dip2*rTd)));
    y1=ti1*(cos((dip1*rTd))*sin((dec1*rTd)));
    y2=ti2*(cos((dip2*rTd))*sin((dec2*rTd)));
    z1=ti1*(sin((dip1*rTd)));
    z2=ti2*(sin((dip2*rTd)));
    h1=ti1*(cos((dip1*rTd)));
    h2=ti2*(cos((dip2*rTd)));

	/*  COMPUTE ANNUAL CHANGE FOR TOTAL INTENSITY  */
    ati = ti2 - ti1;

	/*  COMPUTE ANNUAL CHANGE FOR DIP & DEC  */
    adip = (dip2 - dip1) * 60.;
    adec = (dec2 - dec1) * 60.;

	/*  COMPUTE ANNUAL CHANGE FOR X, Y, Z, AND H */
    ax = x2-x1;
    ay = y2-y1;
    az = z2-z1;
    ah = h2-h1;

    /* deals with geographic and magnetic poles */
      
    if (h1 < 100.0){	/* at magnetic poles */
        dec1 = NaN;
        adec = NaN;
    }
    if (h1 < 1000.0){
        warn_H = 0;
        warn_H_strong = 1;
        warn_H_strong_val = h1;
    }
    else if (h1 < 5000.0 && !warn_H_strong){
          warn_H = 1;
          warn_H_val = h1;
        }
            
     
    /* convert D and I to deg and min */
    if (my_isnan(dec1)) ddeg = dec1; else ddeg=(int)dec1;
    dmin=(dec1-(double)ddeg)*60;
    if (dec1 > 0 && dmin >= 59.5){
        dmin -= 60.0;
		ddeg++;
    }
    if (dec1 < 0 && dmin <= -59.5){
        dmin += 60.0;
        ddeg--;
    }
    if(ddeg!=0) dmin=fabs(dmin);

    if (my_isnan(dip1)) ideg = dip1; else ideg=(int)dip1;
		imin=(dip1-(double)ideg)*60;
		
    if (dip1 > 0 && imin >= 59.5){
          imin -= 60.0;
    ideg++;
    }
   
   if (dip1 < 0 && imin <= -59.5){
        imin += 60.0;
        ideg--;
    }
    
	if(ideg!=0) imin=fabs(imin); 
	  
	inout->F = ti1;
	inout->H = h1;	
	inout->X = x1;	
	inout->Y = y1;	
	inout->Z = z1;	
	inout->decld = ddeg;
	inout->declm = dmin;
	inout->incld = ideg;
	inout->inclm = imin;
}

static void E0000(const char* cof_data, int IENTRY, int *maxdeg, double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv)
{
  static int maxord,i,icomp,n,m,j,D1,D2,D3,D4;
  static double c[13][13],cd[13][13],tc[13][13],dp[13][13],snorm[169],
    sp[13],cp[13],fn[13],fm[13],pp[13],k[13][13],pi,dtr,a,b,re,
    a2,b2,c2,a4,b4,c4,epoch,gnm,hnm,dgnm,dhnm,flnmj,otime,oalt,
    olat,olon,dt,rlon,rlat,srlon,srlat,crlon,crlat,srlat2,
    crlat2,q,q1,q2,ct,st,r2,r,d,ca,sa,aor,ar,br,bt,bp,bpp,
    par,temp1,temp2,parp,bx,by,bz,bh;
  static char model[20], c_new[5];
  const char* c_str;
  static double *p = snorm;
  char answer;
  
  switch(IENTRY){case 0: goto GEOMAG; case 1: goto GEOMG1;}
  
 GEOMAG:
  
/* INITIALIZE CONSTANTS */
  maxord = *maxdeg;
  sp[0] = 0.0;
  cp[0] = *p = pp[0] = 1.0;
  dp[0][0] = 0.0;
  a = 6378.137;
  b = 6356.7523142;
  re = 6371.2;
  a2 = a*a;
  b2 = b*b;
  c2 = a2-b2;
  a4 = a2*a2;
  b4 = b2*b2;
  c4 = a4 - b4;
  
/* READ WORLD MAGNETIC MODEL SPHERICAL HARMONIC COEFFICIENTS */
  c[0][0] = 0.0;
  cd[0][0] = 0.0;
  
  c_str = cof_data;
  cof_data = strchr(cof_data, '\n');
  if (sscanf(c_str,"%lf%s",&epoch,model) < 2) 
   {
       fprintf(stderr, "Invalid header in model data.\n");
       exit(1);
   }

 S3:
  c_str = cof_data;
  cof_data = strchr(cof_data, '\n');
  if (cof_data)
    cof_data++; // skip \n

  int len = cof_data-c_str;

  if (len == 0)
    return;

  if (!c_str)
   goto S4;

/* CHECK FOR LAST LINE IN FILE */
  for (i=0; i<4 && (c_str[i] != '\0'); i++)
    {
      c_new[i] = c_str[i];
      c_new[i+1] = '\0';
    }
  icomp = strcmp("9999", c_new);
  if (icomp == 0) goto S4;
/* END OF FILE NOT ENCOUNTERED, GET VALUES */
  sscanf(c_str,"%d%d%lf%lf%lf%lf",&n,&m,&gnm,&hnm,&dgnm,&dhnm);

  if (n > maxord) goto S4;
  if (m > n || m < 0.0) 
    {
      fprintf(stderr, "Corrupt record in model data.\n");
      exit(1);
    }

  if (m <= n)
    {
      c[m][n] = gnm;
      cd[m][n] = dgnm;
      if (m != 0)
        {
          c[n][m-1] = hnm;
          cd[n][m-1] = dhnm;
        }
    }
  goto S3;

/* CONVERT SCHMIDT NORMALIZED GAUSS COEFFICIENTS TO UNNORMALIZED */
 S4:
  *snorm = 1.0;
  fm[0] = 0.0;
  for (n=1; n<=maxord; n++)
    {
      *(snorm+n) = *(snorm+n-1)*(double)(2*n-1)/(double)n;
      j = 2;
      for (m=0,D1=1,D2=(n-m+D1)/D1; D2>0; D2--,m+=D1)
        {
          k[m][n] = (double)(((n-1)*(n-1))-(m*m))/(double)((2*n-1)*(2*n-3));
          if (m > 0)
            {
              flnmj = (double)((n-m+1)*j)/(double)(n+m);
              *(snorm+n+m*13) = *(snorm+n+(m-1)*13)*sqrt(flnmj);
              j = 1;
              c[n][m-1] = *(snorm+n+m*13)*c[n][m-1];
              cd[n][m-1] = *(snorm+n+m*13)*cd[n][m-1];
            }
          c[m][n] = *(snorm+n+m*13)*c[m][n];
          cd[m][n] = *(snorm+n+m*13)*cd[m][n];
        }
      fn[n] = (double)(n+1);
      fm[n] = (double)n;
    }
  k[1][1] = 0.0;
  
  otime = oalt = olat = olon = -1000.0;
  return;
  
/*************************************************************************/

	GEOMG1:
  
	dt = time - epoch;
	if (otime < 0.0 && (dt < 0.0 || dt > 5.0)){      
		printf("\n\n WARNING - TIME EXTENDS BEYOND MODEL 5-YEAR LIFE SPAN");
    }

  pi = 3.14159265359;
  dtr = pi/180.0;
  rlon = glon*dtr;
  rlat = glat*dtr;
  srlon = sin(rlon);
  srlat = sin(rlat);
  crlon = cos(rlon);
  crlat = cos(rlat);
  srlat2 = srlat*srlat;
  crlat2 = crlat*crlat;
  sp[1] = srlon;
  cp[1] = crlon;

/* CONVERT FROM GEODETIC COORDS. TO SPHERICAL COORDS. */
  if (alt != oalt || glat != olat)
    {
      q = sqrt(a2-c2*srlat2);
      q1 = alt*q;
      q2 = ((q1+a2)/(q1+b2))*((q1+a2)/(q1+b2));
      ct = srlat/sqrt(q2*crlat2+srlat2);
      st = sqrt(1.0-(ct*ct));
      r2 = (alt*alt)+2.0*q1+(a4-c4*srlat2)/(q*q);
      r = sqrt(r2);
      d = sqrt(a2*crlat2+b2*srlat2);
      ca = (alt+d)/r;
      sa = c2*crlat*srlat/(r*d);
    }
  if (glon != olon)
    {
      for (m=2; m<=maxord; m++)
        {
          sp[m] = sp[1]*cp[m-1]+cp[1]*sp[m-1];
          cp[m] = cp[1]*cp[m-1]-sp[1]*sp[m-1];
        }
    }
  aor = re/r;
  ar = aor*aor;
  br = bt = bp = bpp = 0.0;
  for (n=1; n<=maxord; n++)
    {
      ar = ar*aor;
      for (m=0,D3=1,D4=(n+m+D3)/D3; D4>0; D4--,m+=D3)
        {
/*
   COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
   AND DERIVATIVES VIA RECURSION RELATIONS
*/
          if (alt != oalt || glat != olat)
            {
              if (n == m)
                {
                  *(p+n+m*13) = st**(p+n-1+(m-1)*13);
                  dp[m][n] = st*dp[m-1][n-1]+ct**(p+n-1+(m-1)*13);
                  goto S50;
                }
              if (n == 1 && m == 0)
                {
                  *(p+n+m*13) = ct**(p+n-1+m*13);
                  dp[m][n] = ct*dp[m][n-1]-st**(p+n-1+m*13);
                  goto S50;
                }
              if (n > 1 && n != m)
                {
                  if (m > n-2) *(p+n-2+m*13) = 0.0;
                  if (m > n-2) dp[m][n-2] = 0.0;
                  *(p+n+m*13) = ct**(p+n-1+m*13)-k[m][n]**(p+n-2+m*13);
                  dp[m][n] = ct*dp[m][n-1] - st**(p+n-1+m*13)-k[m][n]*dp[m][n-2];
                }
            }
        S50:
/*
    TIME ADJUST THE GAUSS COEFFICIENTS
*/
          if (time != otime)
            {
              tc[m][n] = c[m][n]+dt*cd[m][n];
              if (m != 0) tc[n][m-1] = c[n][m-1]+dt*cd[n][m-1];
            }
/*
    ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS
*/
          par = ar**(p+n+m*13);
          if (m == 0)
            {
              temp1 = tc[m][n]*cp[m];
              temp2 = tc[m][n]*sp[m];
            }
          else
            {
              temp1 = tc[m][n]*cp[m]+tc[n][m-1]*sp[m];
              temp2 = tc[m][n]*sp[m]-tc[n][m-1]*cp[m];
            }
          bt = bt-ar*temp1*dp[m][n];
          bp += (fm[m]*temp2*par);
          br += (fn[n]*temp1*par);
/*
    SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES
*/
          if (st == 0.0 && m == 1)
            {
              if (n == 1) pp[n] = pp[n-1];
              else pp[n] = ct*pp[n-1]-k[m][n]*pp[n-2];
              parp = ar*pp[n];
              bpp += (fm[m]*temp2*parp);
            }
        }
    }
  if (st == 0.0) bp = bpp;
  else bp /= st;
/*
    ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
    GEODETIC COORDINATES
*/
  bx = -bt*ca-br*sa;
  by = bp;
  bz = bt*sa-br*ca;
/*
    COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND
    TOTAL INTENSITY (TI)
*/
  bh = sqrt((bx*bx)+(by*by));
  *ti = sqrt((bh*bh)+(bz*bz));
  *dec = atan2(by,bx)/dtr;
  *dip = atan2(bz,bh)/dtr;
/*
    COMPUTE MAGNETIC GRID VARIATION IF THE CURRENT
    GEODETIC POSITION IS IN THE ARCTIC OR ANTARCTIC
    (I.E. GLAT > +55 DEGREES OR GLAT < -55 DEGREES)

    OTHERWISE, SET MAGNETIC GRID VARIATION TO -999.0
*/
  *gv = -999.0;
  if (fabs(glat) >= 55.)
    {
      if (glat > 0.0 && glon >= 0.0) *gv = *dec-glon;
      if (glat > 0.0 && glon < 0.0) *gv = *dec+fabs(glon);
      if (glat < 0.0 && glon >= 0.0) *gv = *dec+glon;
      if (glat < 0.0 && glon < 0.0) *gv = *dec-fabs(glon);
      if (*gv > +180.0) *gv -= 360.0;
      if (*gv < -180.0) *gv += 360.0;
    }
  otime = time;
  oalt = alt;
  olat = glat;
  olon = glon;
  return;
}

/*************************************************************************/
void geomag(const char* cof_data, int *maxdeg)
{
  E0000(cof_data, 0,maxdeg,0.0,0.0,0.0,0.0,NULL,NULL,NULL,NULL);
}

/*************************************************************************/
void geomg1(const char* cof_data, double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv)
{
  E0000(cof_data, 1,NULL,alt,glat,glon,time,dec,dip,ti,gv);
}
