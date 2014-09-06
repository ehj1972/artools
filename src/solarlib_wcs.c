/* These are WCS-related functions needed for coordinate transformation.
 * 
 * The main list of these can be found in solarlib.h. 
 * 
 * Written by Eric Jones, ejones@sun.stanford.edu
 * 
 * Changelog: 													
 * 
 * 03.22.2014		--		Split off solarlib.c into three files to make 
 * 							organizing things easier. 
 * 
 * 04.24.2014		--		Renamed along_x and family to reflect what they 
 * 							represent a little better, also refined the 
 * 							hcc2sphere function a bit. Added along_theta.	
 * 
 * 06.09.2014		--		Made a number of changes. 
 * 							1. Removed all of the error functions. These will eventually be relocated to solarlib_error.
 * 							2. Added several WCS transformations to get back into pixel space from stonyhurst, etc.
 * 							3. Changed pixel2as to pixel2hpl to denote what it really does.			
 * 
 * */
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "solarlib.h"

/* defines */

#define rsun_ref 696000000.0
#define rad2deg 180.0/M_PI
#define deg2rad M_PI/180.0
#define piov2 M_PI/2.0
#define rad2as rad2deg*3600.0
#define as2rad deg2rad/3600.0

/* Note that the angles theta and phi follow the traditional physics definitions
 * and not the mathematics ones, e.g. here theta is between r and z, and phi is between
 * projection of r onto the x-y plane and the x-axis. */


double along_obs(double theta,double phi,double d_obs){
    /* Finds the LOS weighting of a quantity that is in the observer direction, in 
    heliocentric cartesian coordinates. t,p are in radians (spherical coordinates)
    and d_obs is the keyword dsun_obs. In a conventional spherical coordinate system, this
    would be along the x direction, with z being north and y being "west". */


    /* define variables */
    double num,den;
    double in_obs;
    /* main function */
    num=d_obs-1.0*rsun_ref*sin(theta)*cos(phi);
    den=sqrt(d_obs*d_obs+rsun_ref*rsun_ref-2.0*rsun_ref*d_obs*sin(theta)*cos(phi));
    in_obs=num/den;
    return in_obs;
}

double along_north(double theta,double phi,double d_obs){
    /* Finds the LOS weighting of a quantity that is in the heliographic north direction, in 
    heliocentric cartesian coordinates. t,p are in radians (spherical coordinates)
    and d_obs is the keyword dsun_obs. In a conventional spherical coordinate system, this
    would be along the z direction, with x being out of the page and y being "west". */


    /* define variables */
    double num,den;
    double in_north;
    /* main function */
    num=(-1.0*rsun_ref*cos(theta));
    den=sqrt(d_obs*d_obs+rsun_ref*rsun_ref-2.0*rsun_ref*d_obs*sin(theta)*cos(phi));
    in_north=num/den;
    return in_north;
}

double along_west(double theta,double phi,double d_obs){
    /* Finds the LOS weighting of a quantity that is in the heliographic west direction, in 
    heliocentric cartesian coordinates. t,p are in radians (spherical coordinates)
    and d_obs is the keyword dsun_obs. In a conventional spherical coordinate system, this
    would be along the y direction, with x being out of the page and z being north. */


    /* define variables */
    double num,den;
    double in_west;
    /* main function */
    num=-1.0*rsun_ref*sin(theta)*sin(phi);
    den=sqrt(d_obs*d_obs+rsun_ref*rsun_ref-2.0*rsun_ref*d_obs*sin(theta)*cos(phi));
    in_west=num/den;
    return in_west;
}

double along_theta(double theta,double phi,double d_obs){
    /* Finds the LOS weighting of a quantity that is in the theta direction, in 
    standard spherical coordinates. t,p are in radians (spherical coordinates)
    and d_obs is the keyword dsun_obs */


    /* define variables */
    double num,den;
    double in_theta;
    /* main function */
    num=d_obs*cos(theta)*cos(phi);
    den=sqrt(d_obs*d_obs+rsun_ref*rsun_ref-2.0*rsun_ref*d_obs*sin(theta)*cos(phi));
    in_theta=num/den;
    return in_theta;
}

double along_phi(double theta,double phi,double d_obs){
    /* Finds the LOS weighting of a quantity that is in the phi direction, in 
    standard spherical coordinates. t,p are in radians (spherical coordinates)
    and d_obs is the keyword dsun_obs */


    /* define variables */
    double num,den;
    double in_phi;
    /* main function */
    num=-1.0*d_obs*sin(phi);
    den=sqrt(d_obs*d_obs+rsun_ref*rsun_ref-2.0*rsun_ref*d_obs*sin(theta)*cos(phi));
    in_phi=num/den;
    return in_phi;
}

double along_r(double theta,double phi,double d_obs){
    /* Finds the LOS weighting of a quantity that is in the radial direction, in 
    standard spherical coordinates. t,p are in radians (spherical coordinates)
    and d_obs is the keyword dsun_obs */


    /* define variables */
    double num,den;
    double in_r;
    /* main function */
    num=d_obs*sin(theta)*cos(phi)-rsun_ref;
    den=sqrt(d_obs*d_obs+rsun_ref*rsun_ref-2.0*rsun_ref*d_obs*sin(theta)*cos(phi));
    in_r=num/den;
    return in_r;
}


void pixel2hpl(double x, double y, double crota2, double crpix1, double crpix2, double cdelt1, double cdelt2, double* tx, double* ty){

    /* Converts a set of full-disk pixel coordinates from HMI into pseudo-angle HPLN-TAN/HPLT-TAN
    coordinates, relative to the FITS center-of-image crpix1/2 location and rotating 
    counterclockwise as per the crota2 keyword.
    
    NOTE: Assumes HPLN-TAN/HPLT-TAN projection.	*/	
	
	double x0,y0;
	double trad;
	double nx,ny;
	
	x0=crpix1-1;
	y0=crpix2-1;
	trad=deg2rad*crota2;
	
	rotate(x-x0,y-y0,trad,&nx,&ny);
	*tx=cdelt1*nx;
	*ty=cdelt2*ny;
	return;
}

void hpl2pixel(double tx, double ty, double crota2, double crpix1, double crpix2, double cdelt1, double cdelt2, double *x, double *y){
	
	/* Converts pseudo-angle HPLN-TAN/HPLT-TAN into full-disk coordinates, de-rotating as we go. */
	
	double x0,y0;
	double trad;
	double nx,ny;
	double xx,yy;
	
	trad=deg2rad*crota2;
	x0=crpix1-1;
	y0=crpix2-1;
	
	nx=tx/cdelt1;
	ny=ty/cdelt2;
	rotate(nx,ny,-1.0*trad,&xx,&yy);
	*x=xx+x0;
	*y=yy+y0;
	
	return;
	
}
	
void hpl2hpc(double x, double y, double* tx, double* ty){
	 
	/* Converts HPLN/HPLT-TAN pseudo angle coordinates to helioprojective-cartesian
    coordinates.	*/
    
    double xx,yy;
    double g;
    
    xx=as2rad*x;
    yy=as2rad*y;
    g=yy/sqrt(1.0+xx*xx);
    *tx=rad2as*atan(xx);
    *ty=rad2as*atan(g);
	return;
}

void hpc2hpl(double tx, double ty, double *x, double *y){

	/* Converts helioprojective cartesian coordinates to HPLN/HPLT-TAN ones */
	
	double ttx,tty;
	
	ttx=as2rad*tx;
	tty=as2rad*ty;
	
	*x=rad2as*tan(ttx);
	*y=rad2as*tan(tty)/cos(ttx);
	
	return;
}

void hpc2hcc(double ttx, double tty, double dsun_obs, double* x, double* y, double* z){
	
	/* Converts helioprojective-cartesian coordinates into heliocentric-cartesian 
    coordinates. Returns -1000's for x,y,z if there is no solution.	*/
    
    double tx,ty;
    double A,B,C;
    double AA,BB,d;
    
    tx=as2rad*ttx;
    ty=as2rad*tty;
    
    A=1.0;
    B=-2.0*dsun_obs*cos(tx)*cos(ty);
    C=pow(dsun_obs,2.0)-pow(rsun_ref,2.0);
    AA=4.0*A*C;
    BB=B*B;
    
    /* Make sure things don't get imaginary */
    
    if (AA > BB) {
		
		*x=*y=*z=-1000;
	
	} else {
		
		/* The negative solution for the quadratic gives the correct physics */
		
		d=(-1.0*B-sqrt(BB-AA))/2.0;
			
		*x=d*cos(ty)*sin(tx);
		*y=d*sin(ty);
		*z=dsun_obs-d*cos(ty)*cos(tx);
	
	}
	
	return;

}

void hcc2hpc(double x, double y, double z, double dsun_obs, double *tx, double *ty){
	
	/* Converts heliocentric cartesian to helioprojective cartesian */
	
	double d;
	
	d=sqrt(x*x+y*y+pow(dsun_obs-z,2.0));
	*tx=rad2as*arg(dsun_obs-z,x);
	*ty=rad2as*asin(y/d);
	
	return;
	
}

void hcc2stony(double x, double y, double z, double P, double B, double* lat, double* lon){
	
	/*Converts heliocentric cartesian coordinates to Stonyhurst heliographic 
    coordinates.	*/
    
    double P0,B0;
    double r;
    
    P0=deg2rad*P;
    B0=deg2rad*B;
    r=sqrt(x*x+y*y+z*z);
    *lat=rad2deg*(asin((y*cos(B0)+z*sin(B0))/r));
    *lon=rad2deg*(P0+arg(z*cos(B0)-y*sin(B0),x));
    
	return;	
}

void stony2hcc(double lat, double lon, double P, double B, double *x, double *y, double *z){

	/* Converts Stonyhurst heliographic coordinates into heliocentric cartesian ones */

    double P0,B0;
    double latitude,longitude;
    
    P0=deg2rad*P;
    B0=deg2rad*B;
    latitude=deg2rad*lat;
    longitude=deg2rad*lon;
    
    *x=rsun_ref*cos(latitude)*sin(longitude-P0);
    *y=rsun_ref*(sin(latitude)*cos(B0)-cos(latitude)*cos(longitude-P0)*sin(B0));
    *z=rsun_ref*(sin(latitude)*sin(B0)+cos(latitude)*cos(longitude-P0)*cos(B0));

	return;

}

void hpl2stony(double x, double y, double dsun_obs, double P0, double B0, double* lat, double* lon){
	
	/*Converts HPLN/HPLT-TAN coordinates to Stonyhurst coordinates. Values of dsun_obs, 
    x, and y that are off-disk will cause an error and will return as lat,lon = -1000,-1000.	*/
    
    double tx,ty;
    double true_x,true_y,true_z;
    
    hpl2hpc(x,y,&tx,&ty);
    hpc2hcc(tx,ty,dsun_obs,&true_x,&true_y,&true_z);
    
    if ((true_x == -1000)||(true_y==-1000)||(true_z==-1000)){
	
		*lat=*lon=-1000;
	
	} else {
	
		hcc2stony(true_x,true_y,true_z,P0,B0,lat,lon);
	
	}
	
	return;
}

void stony2hpl(double lat, double lon, double dsun_obs, double P0, double B0, double *x, double *y){
	
	/* Converts coordinates in Stonyhurst to HPLN/HPLT-TAN coordinates. */
	
	double tx,ty;
	double true_x,true_y,true_z;
	
	if ((lat == -1000)||(lon == -1000)){
		
		*x=-10000;
		*y=-10000;
	
	} else {
		
		stony2hcc(lat,lon,P0,B0,&true_x,&true_y,&true_z);
		hcc2hpc(true_x,true_y,true_z,dsun_obs,&tx,&ty);
		hpc2hpl(tx,ty,x,y);
		
	}
	return;
}


void hcc2sphere(double helio_west,double helio_north,double observer,double *theta,double *phi){
    /* This converts heliocentric cartesian to spherical polar coordinates,
     * such that theta is the angle between the r vector and solar north,
     * and phi is between the observer direction and the projection of r 
     * onto the observer-solar west plane. Think standard spherical polar,
     * but adapted to Thompson's 2005 paper.	*/

    /* define some variables */
    double f,g;
    double x,y,z;
    /* conversions from thompson to standard spherical polar */
    
    z=helio_north;
    y=helio_west;
    x=observer;
    
    g=z/rsun_ref;
    f=fabs(y/x);	/*fabs for quadrant stuff */
 
    /* Handle quadrant issues with tan */
    
    if ((x > 0)&&(y > 0)){
		/* first quad */
		*phi=atan(f);
	}
	
    if ((x > 0)&&(y < 0)){
		/* fourth quad */
		*phi=2.0*M_PI-atan(f);
	}
	
    if ((x < 0)&&(y > 0)){
		/* second quad */
		*phi=M_PI-atan(f);
	}
	
    if ((x < 0)&&(y < 0)){
		/* third quad */
		*phi=M_PI+atan(f);
	}
	
	/* Special cases */
	
	if ((x == 0)&&(y == 0)){
		*phi=0;
	}
	
	if ((x > 0)&&(y == 0)){
		*phi=0;
	}
	
	if ((x == 0)&&(y > 0)){
		*phi=M_PI/2.0;
	}
	
	if ((x < 0)&&(y == 0)){
		*phi=M_PI;
	}
	
	if ((x == 0)&&(y < 0)){
		*phi=(3.0/2.0)*M_PI;
	}
	
	*theta=acos(g);
	
    return;
}

void CEA2latlon(double cea_lat, double cea_lon, double lon_cm, double *lat_deg, double *lon_deg){
	
	/* Converts coordinates in Lambert CEA to standard latitude and longitude */
	
	*lon_deg=cea_lon + lon_cm;
	*lat_deg=rad2deg*asin(cea_lat);
	
	return;
	
}

	 

    
    
    
	




	
	
	
	

