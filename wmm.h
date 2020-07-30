/***********************************************************/
/* All rights reserved by Amir Hossein Alikhah Mishamandani*/
/******* Created for ADCS system of SSS-P1 Satellite *******/
/*********************** 2018 - 2019 ***********************/
/************ World Magnetic Model 2015 - 2020 *************/
/*********** Original Program By: Dr. John Quinn ***********/
/********** Corrected by: Stefan Maus & Manoj Nair *********/
/***********************************************************/

#ifndef __WMM_H__
#define __WMM_H__
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    double X;		//nT		
    double Y;		//nT		
    double Z;		//nT		
    double decld;	//Deg
	double declm;	//Min
    double incld;	//Deg
	double inclm;	//Min
    double H;		//nT
    double F;		//nT
} geomag_vector;

void WorldMagneticModel(const char* cof_data, geomag_vector *inout, double dlat,double dlon,double altm,double time);
//dlat is LATITUDE (IN DECIMAL DEGREES): North latitude positive & South latitude negative. (i.e. 25.5 for 25 degrees 30 minutes north.)
//dlon is LONGITUDE (IN DECIMAL DEGREES): East longitude positive & West negative. (i.e.- 100.0 for 100.0 degrees west.)
//altm is ALTITUDE (IN KILOMETERS): ABOVE WGS84 ELLIPSOID.
//time (IN DECIMAL YEAR): 2015 - 2020


#ifdef __cplusplus
}
#endif

#endif