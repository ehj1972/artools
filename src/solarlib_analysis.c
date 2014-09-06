/* These are functions specifically used during the data analysis sections
 * of code.
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
 * 03.30.2014		--		Added connected_component_detect, which finds
 * 							discrete regions within an image and labels
 * 							them, returning the number of discrete elements
 * 							within the image.
 * 				
 * 06.15.2014		--		General update. Since last, a number of things have
 * 							been added. This time, adding pixel_area.
 * 
 * 08.09.2014		--		Changed spotbounds so that cutoffs are function arguments
 * 							instead of being hard-coded in.
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

double limb_darken(double mu){
	
	/* Returns the continuum intensity at mu of a HMI-observed quiet sun, normalized
    to a 10x10 pixel patch centered at crpix1,crpix2. Fit coefficients were determined using 
    the polynomial form in Hestroffer and Magnan, 1998. */
	
	double a0,a1,a2;		/* Coefficients from fit described above 	*/
	double I;				/* Normalized quiet Sun intensity			*/
	a0=0.36795341;
	a1=0.8794135;
	a2=-0.23355535;
	I=a0+a1*mu+a2*mu*mu;
	return I;
}

void spotbounds(double** cdat, int naxis1, int naxis2, double* umb, double* pen, double umb_sig, double pen_sig){
	
	/* Calculates the cutoff levels for the umbra and the penumbra within a 
    sunspot. Sunspot is assumed to be umb_sig away from background, with
    umbra being pen_sig away from the average of the sunspot intensity. 0 values are 
    assumed to be bad and not included in mean or stdev calculations.
    (Verbeeck et al, 2011) */

    double mean, stdev;
    double smean,sstdev;
    int i,j;
    double counter;
    
    counter=0.0;
    mean=smean=0.0;
    stdev=sstdev=0.0;
    
    *umb=0.0;
    *pen=0.0;
    
    /* First calculate the mean and stdev of the patch */
    
    for (i=0;i<naxis1;i++){
		for (j=0;j<naxis2;j++){
			if (cdat[j][i] > 0){
				mean+=cdat[j][i];
				counter++;
			}
		}
	}
	
	mean/=counter;
	
    for (i=0;i<naxis1;i++){
		for (j=0;j<naxis2;j++){
			if (cdat[j][i] > 0){
				stdev+=pow(cdat[j][i]-mean,2.0);
			}
		}
	}
	
	stdev=sqrt(stdev/(counter-1.0));
	
	/* Now find the umbra and penumbra cutoffs */
	
	*pen=mean-umb_sig*stdev;
	counter=0.0;
	
    for (i=0;i<naxis1;i++){
		for (j=0;j<naxis2;j++){
			if ((cdat[j][i] > 0)&&(cdat[j][i] <= *pen)){
				smean+=cdat[j][i];
				counter++;
			}
		}
	}
	
	smean/=counter;
	
    for (i=0;i<naxis1;i++){
		for (j=0;j<naxis2;j++){
			if ((cdat[j][i] > 0)&&(cdat[j][i] <= *pen)){
				sstdev+=pow(cdat[j][i]-smean,2.0);
			}
		}
	}
	
	sstdev=sqrt(sstdev/(counter-1.0));
	*umb=smean-pen_sig*sstdev;
	
	return;
}

void arclength(double lon_1,double lon_0,double lat_1,double lat_0, double* arclen){
	
	/* Calculates the arc length between two points on the surface of the Sun.
    Assumes that the radius of the Sun is constant. Results are returned in meters, 
    and centroid pairs that contain a -1000 return a 0 to show that either the 
    sunspots are below threshold or one of the polarities is.	*/
    
    if ((lat_1 == -1000)||(lat_0 == -1000)||(lon_1 == -1000)||(lon_0 == -1000)){
		
		*arclen=0.0;
	
	} else {
		
		*arclen=rsun_ref*acos(cos(deg2rad*fabs(lat_1-lat_0))*cos(deg2rad*fabs(lon_1-lon_0)));

	}
	
	return;
}

void arclength_e(double lon_1,double lon_0,double lat_1,double lat_0, double dlon_1, double dlon_0, double dlat_1, double dlat_0, double* arclen, double *darclen){
	
	/* Calculates the arc length and error between two points on the surface of the Sun.
    Assumes that the radius of the Sun is constant. Results are returned in meters, 
    and centroid pairs that contain a -1000 return a 0 to show that either the 
    sunspots are below threshold or one of the polarities is.	*/
    
    double g;
    double dgdlon_0,dgdlon_1;
    double dgdlat_0,dgdlat_1;
    double darcdlon_0,darcdlon_1;
    double darcdlat_0,darcdlat_1;
    
    if ((lat_1 == -1000)||(lat_0 == -1000)||(lon_1 == -1000)||(lon_0 == -1000)){
		
		*arclen=0.0;
		*darclen=0.0;
	
	} else {
		
		*arclen=rsun_ref*acos(cos(deg2rad*fabs(lat_1-lat_0))*cos(deg2rad*fabs(lon_1-lon_0)));
            
		g=cos(deg2rad*(lat_1-lat_0))*cos(deg2rad*(lon_1-lon_0));
		
		dgdlon_0=cos(deg2rad*(lat_1-lat_0))*sin(deg2rad*(lon_1-lon_0));
		dgdlon_1=-1.0*cos(deg2rad*(lat_1-lat_0))*sin(deg2rad*(lon_1-lon_0));
		dgdlat_0=sin(deg2rad*(lat_1-lat_0))*cos(deg2rad*(lon_1-lon_0));
		dgdlat_1=-1.0*sin(deg2rad*(lat_1-lat_0))*cos(deg2rad*(lon_1-lon_0));

		darcdlon_0=(-1.0*rsun_ref/sqrt(1.0-g*g))*dgdlon_0;
		darcdlon_1=(-1.0*rsun_ref/sqrt(1.0-g*g))*dgdlon_1;
		darcdlat_0=(-1.0*rsun_ref/sqrt(1.0-g*g))*dgdlat_0;
		darcdlat_1=(-1.0*rsun_ref/sqrt(1.0-g*g))*dgdlat_1;
		
		*darclen=sqrt(pow(darcdlon_0*dlon_0,2.0) + pow(darcdlon_1*dlon_1,2.0) + pow(darcdlat_0*dlat_0,2.0) + pow(darcdlat_1*dlat_1,2.0));		

	}
	
	return;
}

void tilt_angle(double lon_1,double lon_0,double lat_1,double lat_0, double* tilt){
	
	/* Calculates the tilt angle of a sunspot group, based upon the latitudes and longitudes of the
    polarity centroids. Returns 0 if one or more of the coordinates is bad (defined as -1000). */ 
    
    double g;
    
    if ((lat_1 == -1000)||(lat_0 == -1000)||(lon_1 == -1000)||(lon_0 == -1000)){
		
		*tilt=0.0;
		
	} else {
		
		/* Sun rotates from small angle to large angle, so we can use that to determine sign of tilt. */
		
		g=(double)fabs((lat_1-lat_0)/(lon_1-lon_0));
		
		*tilt=rad2deg*atan(g);
		
		/* Ok. Some checking here. If lat > 0 then it's a southern hemisphere location, otherwise northern */
		
		if ((lat_0 > 0)&&(lat_1 > 0)){
			
			/* Northern Hemisphere. Here, bigger latitude means bigger tilt. */		
		
			if (lon_1 > lon_0){
			
				/* lon_1 leads, therefore if lat_1 > lat_0 it's tilting contrary to Joy's law */
			
				if (lat_1 > lat_0){
				
					/* Negative tilt means leading centroid lower than following, does not obey Joy's Law.  */
				
					*tilt*=-1.0;
				
				}
			
			} else {
			
				/* lon_0 leads. Same argument as above */
			
				if (lat_0 > lat_1){
				
					/* Negative tilt means leading centroid lower than following, does not obey Joy's Law.  */
				
					*tilt*=-1.0;
				
				}
			
			}
			
		}
		
		if ((lat_0 < 0)&&(lat_1 < 0)){
			
			/* Southern Hemisphere. Here, smaller (more negative) latitude means bigger tilt. */		
		
			if (lon_1 > lon_0){
			
				/* lon_1 leads, therefore if lat_0 > lat_1 it's tilting contrary to Joy's law */
			
				if (lat_1 < lat_0){
				
					/* Negative tilt means leading centroid lower than following, does not obey Joy's Law.  */
				
					*tilt*=-1.0;
				
				}
			
			} else {
			
				/* lon_0 leads. Same argument as above */
			
				if (lat_0 < lat_1){
				
					/* Negative tilt means leading centroid lower than following, does not obey Joy's Law.  */
				
					*tilt*=-1.0;
				
				}
			
			}	
				
		}

	}
	
	return;
	
}

void tilt_angle_e(double lon_1,double lon_0,double lat_1,double lat_0, double dlon_1, double dlon_0, double dlat_1, double dlat_0, double* tilt, double* dtilt){
	
	/* Calculates the tilt angle of a sunspot group and error, based upon the latitudes and longitudes of the
    polarity centroids. Returns 0 if one or more of the coordinates is bad (defined as -1000). */ 
    
    double g;
    double dgdlat_1,dgdlat_0;
    double dgdlon_1,dgdlon_0;
    double dtiltdlat_1,dtiltdlat_0;
    double dtiltdlon_1,dtiltdlon_0;
    
    if ((lat_1 == -1000)||(lat_0 == -1000)||(lon_1 == -1000)||(lon_0 == -1000)){
		
		*tilt=0.0;
		*dtilt=0.0;
		
	} else {
		
		/* Sun rotates from small angle to large angle, so we can use that to determine sign of tilt. */
		
		g=(double)fabs((lat_1-lat_0)/(lon_1-lon_0));
		
		*tilt=rad2deg*atan(g);

		/* Error */
		
		dgdlat_1=1.0/(lon_1-lon_0);
		dgdlat_0=-1.0/(lat_1-lat_0);
		dgdlon_1=-1.0*g/(lon_1-lon_0);
		dgdlon_0=g/(lat_1-lat_0);

		dtiltdlat_1=(1.0/(1.0+g*g))*dgdlat_1;
		dtiltdlat_0=(1.0/(1.0+g*g))*dgdlat_0;
		dtiltdlon_1=(1.0/(1.0+g*g))*dgdlon_1;
		dtiltdlon_0=(1.0/(1.0+g*g))*dgdlon_0;
 
		*dtilt=rad2deg*sqrt(pow(dtiltdlat_1*dlat_1,2.0) + pow(dtiltdlat_0*dlat_0,2.0) + pow(dtiltdlon_1*dlon_1,2.0) + pow(dtiltdlon_0*dlon_0,2.0));		
		
		if (lon_1 > lon_0){
			
			if (lat_1 < lat_0){
				
				/* Negative tilt means leading centroid lower than following, does not obey Joy's Law.  */
				
				*tilt*=-1.0;
				
			}
			
		} else {
			
			if (lat_0 < lat_1){
				
				/* Negative tilt means leading centroid lower than following, does not obey Joy's Law.  */
				
				*tilt*=-1.0;
				
			}	
			
		}	
				
	}

	return;
	
}

double centroid_error(double* x, double* m_err, int length, double Xc, double M){

	/*This is an error calculation routine for determining the standard error on centroid calculations
	  of the form:
        
      X = sum(m_i*x_i)/sum(m_i)
        
      where the error is standard and given for m_i, which is analagous to mass but could be basically
      any quantity, and x_i is position. The error propagation equation for this using differential 
      calculus is:
        
      err(X)^2 = sum((((1/M)*(x_i-X))^2)*err(m_i)^2)
    
      where M is sum(m_i). In the event x has a length of 0, -1 is returned as the error. */
      
	double x_err;
	int i;
      
	x_err=0.0;
      
	for (i=0;i<length;i++){
		  
		x_err+=pow((x[i]-Xc)/M,2.0)*pow(m_err[i],2.0);
		
	}
	
	x_err=sqrt(x_err);
	
	return x_err;
	
}

double br_error(double c_min, double c_quiet, double c_min_err, double c_quiet_err){
	
	/*This is the error calculation routine for determining the standard error on the brightness ratio, 
      which has the form:

      BR=c_min/c_quiet

      where the error is standard and given for c_min and c_quiet. The error propagation equation for this 
      using differential calculus is:

      err(BR)^2 = ((1/c_quiet)*c_quiet_err)^2 + ((-BR/c_quiet)*c_min_err)^2 */
      
	double br,br_err;
	
	br=c_min/c_quiet;
	
	br_err=sqrt(pow(c_quiet_err/c_quiet,2.0) + pow(br*c_min_err/c_quiet,2.0));
	
	return br_err;

}

double v_t(double latitude){
    /* Calculates the tangential velocity of the photosphere at latitude. 
    Coefficients in are from Snodgrass 1984a (see Beck's 1999 summary), using 
    photospheric spectral data. */

    /* define variables */
    double lat,theta;
    double A,B,C;
    double v_tan;
    /* main function */
    A=2.851;
    B=-0.343;
    C=-0.474;
    lat=latitude*deg2rad;
    theta=piov2-lat;
    v_tan=omega(lat,A,B,C)*rsun_ref*sin(theta);
    return v_tan;
}
      
double is_normal(double t, double p, double d_obs){
	
	/* Calculates the scaling needed to turn a LOS quantity into the true
	 * value assuming that the quantity is in fact normal to the solar surface */
	
	
	double num,den;
	double cos_g;
	
	num=(d_obs*sin(t)*cos(p))-rsun_ref;
    den=sqrt(pow(d_obs,2.0)+pow(rsun_ref,2.0)-2.0*rsun_ref*d_obs*sin(t)*cos(p));
    cos_g=num/den;
    return 1.0/cos_g;
}

int connected_component_detect(int **data, int naxis1, int naxis2){
	
	/* Calculates the number of discrete regions within an image. The 
	 * image passed must be binary, with background set to 0 and foreground
	 * set to 1. */
	
	int **labeling;
	int i,j;
	int newlabel;
	double fg_value;
	int temp_label;


	/* Create and zero out the label array */
	
	labeling=malloc(naxis2*sizeof(int *));
	
	for (i=0;i<naxis2;i++){
		labeling[i]=malloc(naxis1*sizeof(int));
	}	
	
	for (i=0;i<naxis2;i++){
		for (j=0;j<naxis1;j++){
			labeling[i][j]=0;
		}
	}	

	newlabel=1;
	fg_value=1.0;

	/* First iteration, set up the labels */
	
	for (i=0;i<naxis2;i++){
		
		for (j=0;j<naxis1;j++){
			
			/* Check to see if the pixel is background or foreground */
			
			if (data[i][j] == fg_value){
					
				/* Do any neighbors have a label? If so, find the smallest nonzero 
				 * label. If not, assign one. */
					
				temp_label=check_labels(i,j,labeling,naxis2,naxis1);
				
				/* This returns either a number or a 0. If it's 0, we need to add a new 
				 * label to this pixel. If it's a number, then we should set the label 
				 * to this value. */
				 
				if (temp_label > 0){
					
					labeling[i][j]=temp_label;
				
				} else {
					
					labeling[i][j]=newlabel;
					newlabel++;
				
				}
					
			}
		}
	}
	
	/* Second iteration, fix any problem labels */
	
	for (i=0;i<naxis2;i++){
		
		for (j=0;j<naxis1;j++){
			
			/* Check to see if the pixel is background or foreground */
			
			if (data[i][j] == fg_value){
				
				/* Do any neighbors have a label? If so, find the smallest nonzero 
				 * label. If not, assign one. */
					
				temp_label=check_labels(i,j,labeling,naxis2,naxis1);
				
				/* This should now return a nonzero number, which we assign to the 
				 * label of the pixel. */
				 
				labeling[i][j]=temp_label;
				
			}
		}
	}
	
	free(labeling);
	
	return (newlabel-1);
}

double pixel_area(double lat, double lon, double dpix){
	
	/* Determines the area of a pixel on the Sun, centered on lat and lon
	 * with an angular size on a side of dpix. */
	 
	double area;
	double ti,tf,pi,pf;
	double t,p;
	double da;
	
	da=deg2rad*dpix/2.0;
	
	/* Standard spherical coordinates - only positive angles */
	
	if (lon < 0){
		
		p=deg2rad*(360+lon);
		
	} else {
		
		p=deg2rad*lon;
	}
	
	/* Need to go from latitude where equator=0 to latitude where north pole=0 */
	
	t=deg2rad*(90.0-lat);
	
	/* Now the area */
	
	pi=p-da;
	pf=p+da;
	ti=t-da;
	tf=t+da;
		
	area=-1.0*(pow(rsun_ref,2.0))*(pf-pi)*(cos(tf)-cos(ti));
	
	return area;
	
}
