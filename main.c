/***********************************************************/
/* All rights reserved by Amir Hossein Alikhah Mishamandani*/
/******* Created for ADCS system of SSS-P1 Satellite *******/
/*********************** 2018 - 2019 ***********************/
/******************* Main ADCS Platform ********************/
/***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "wmm.h"

int main(){
	
	geomag_vector now;
	time_t my_time;
	struct tm * timeinfo; 
	time (&my_time);
	timeinfo = localtime (&my_time);
	double lat = 15.8;
	double lon = 15.6;

	char buffer[1024*1024];
	FILE *wmmtemp = fopen("wmm.cof","rb");
	if (wmmtemp == NULL){
		fprintf(stderr, "Error opening model file WMM.COF\n");
		exit(1);
	}

	fread(buffer, 1, sizeof(buffer), wmmtemp);
	fclose(wmmtemp);


	WorldMagneticModel(buffer, &now, lat, lon, 500, timeinfo->tm_year+1900);// Lat, Lon, Alt, Time
	printf("\nInitializing the World Magnetic Model Test:");
	printf("\nThe latitude = %f", lat);
	printf("\nThe longitude = %f", lon);
	printf("\nThe Altitude = 500");
	printf("\nThe year = %d",timeinfo->tm_year+1900);
	printf("\nThe World Magnetic Model Test successfull ...\n");
	printf("\nX - North Component of the geomagnetic field:	%lf	nT", now.X);
	printf("\nY - East Component of the geomagnetic field:	%lf	nT", now.Y);
	printf("\nZ - Vertical Component of the geomagnetic field:	%lf	nT", now.Z);
	printf("\nF - Total Intensity of the geomagnetic field:	%lf	nT", now.F);
	printf("\nH - Horizontal Intensity of the geomagnetic field:	%lf	nT", now.H);
	printf("\nI (DIP) - Geomagnetic Inclination:	%lf	deg	%lf	min", now.incld, now.inclm);
	printf("\nD (DEC) - Geomagnetic Declination (Magnetic Variation):	%lf	deg	%lf	min\n", now.decld, now.declm);
}
