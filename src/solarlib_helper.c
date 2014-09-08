/* These are helper functions that don't fall into the wcs or analysis category.
 * The main list of these can be found in solarlib.h. 
 * 
 * Written by Eric Jones, ejones@sun.stanford.edu
 * 
 * Changelog: 													
 * 
 * 03.22.2014		--		Split off solarlib.c into three files to make 
 * 							organizing things easier. 		
 * 
 * 03.30.2014		--		Added check_labels, needed by connected_component_detect
 * 							in the analysis section.
 * 	
 * 05.29.2014		--		Added int_to_binary to handle checking error codes in info_map and elsewhere.
 * 
 * 06.02.2014		--		Added bilinear_interpolate and supporting functions.
 * 
 * 09.08.2014		--		Moved omega function from this file to solarlib_analysis.
 * 
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

void rotate(double x, double y, double theta, double* rx, double* ry){
	
	/* Rotates a vector (or set of vectors) with lengths x,y by angle theta (or sets of angles 
    corresponding to the vector set), giving new lengths. Theta is defined contrary to the standard 
    physics definition, with CW being positive rotation, meeting HMI's definition.	*/
    
    *rx=x*cos(theta)-y*sin(theta);
    *ry=x*sin(theta)+y*cos(theta);
    return;
}

double arg(double x, double y){
	
	/* This is equivalent to the trig "arg" function, as defined by footnotes in
    Thompson, 2005 for solar coordinate systems. It restricts things between -180 and
    180 degrees to get rid of coordinate ambiguities. Treat as arctan(y/x).	*/
    
    double angle;
    
    /* Quadrant handling */
    
    if ((x>0)&&(y>0)){
		angle=atan(y/x);
	}
	
    if ((x<0)&&(y>0)){
		angle=M_PI-atan(fabs(y/x));
	}
	
    if ((x>0)&&(y<0)){
		angle=-1.0*atan(fabs(y/x));
	}
	
    if ((x<0)&&(y<0)){
		angle=atan(fabs(y/x))-M_PI;
	}
	
	/* Special cases */
	
    if ((x==0)&&(y>0)){
		angle=M_PI/2.0;
	}
    
    if ((x==0)&&(y<0)){
		angle=3.0*M_PI/2.0;
	}
    
    if ((x>0)&&(y==0)){
		angle=0.0;
	}
	
	if ((x<0)&&(y==0)){
		angle=M_PI;
	}
	
	if ((x==0)&&(y==0)){
		angle=0;
	}
	return angle;
}

long isotime(char* t_obs){

	/* Takes a 32-character T_OBS string and converts it to seconds, referenced to the UNIX epoch. */
	
	int year,month,day;
	int hour,minute,second;
	char s_year[5],s_month[3],s_day[3];
	char s_hour[3],s_minute[3],s_second[3];
	
	struct tm unix_time;
	time_t t_obs_secs;
	
	/* Copy from t_obs */
	
	strncpy(s_year,t_obs,4);
	strncpy(s_month,t_obs+5,2);
	strncpy(s_day,t_obs+8,2);
	strncpy(s_hour,t_obs+11,2);
	strncpy(s_minute,t_obs+14,2);
	strncpy(s_second,t_obs+17,2);
	
	/* Handle EOL garbage */
	
	s_year[4]='\0';
	s_month[2]='\0';
	s_day[2]='\0';
	s_hour[2]='\0';
	s_minute[2]='\0';
	s_second[2]='\0';
	
	/* Convert to integers */
	
	sscanf(s_year,"%d",&year);
	sscanf(s_month,"%d",&month);
	sscanf(s_day,"%d",&day);
	sscanf(s_hour,"%d",&hour);
	sscanf(s_minute,"%d",&minute);
	sscanf(s_second,"%d",&second);
	
	/* Place in time structure, find seconds since epoch */
	
	unix_time.tm_year=year-1900;
	unix_time.tm_mon=month-1;
	unix_time.tm_mday=day;
	unix_time.tm_hour=hour;
	unix_time.tm_min=minute;
	unix_time.tm_sec=second;
	unix_time.tm_isdst=0;
	
	t_obs_secs=mktime(&unix_time);
	
	return (long)t_obs_secs;
	
}

int check_labels(int i, int j, int** label_array, int rows, int cols){
	
	/* Checks the labels of nearest pixels and returns the lowest nonzero label. 
	 * Returns 0 if no neighbors are labeled, including the original pixel. */
	
	int k,l;
	int x,y;
	int min_value;
	int current_label;
	
	current_label=label_array[i][j];
	min_value=10000;
	
	for (k=-1;k<=1;k++){
		for (l=-1;l<=1;l++){
			
			/* Sanity check to make sure we stay within the actual data */
			x=i+k;
			y=j+l;
		
			if ((x < rows)&&(x >= 0)&&(y < cols)&&(y >= 0)){
			
				/* Find the smallest nonzero label, including that of the pixel we are looking at */
			
				if ((label_array[x][y] > 0)&&(label_array[x][y] < min_value)){
				
					min_value=label_array[x][y];
			
				}
				
			}
		}	
	}
	
	/* Ok, so here are the conditions. If min_value hasn't changed, all of the
	 * neighboring pixels are 0. If current_label is 0 as well, we need to add
	 * a new label. To show this, we pass back a 0.
	 * 
	 * If min_value has changed, then we need to compare it to current_label and 
	 * choose the lowest nonzero label, which should automatically have been done
	 * in the above if statement.	*/
				
	if ((min_value == 10000)&&(current_label == 0)){
				
		min_value=0;
					
	}
	
	return min_value;
}

void draw_centroid(double x, double y, double radius, int naxis1, int naxis2, int **data, int color){
	
	/* This function draws a grayscale-filled circle of radius "radius", centered at x,y within the
	 * dataset data. This is for jpeg output. */
	 
	 int i,j;
	 double r;
	 
	 for (j=0;j<naxis2;j++){
		 
		 for (i=0;i<naxis1;i++){
			 
			 r=sqrt(pow((double)i-x,2.0)+pow((double)j-y,2.0));
			 
			 if (r <= radius){
				 
				 data[j][i]=color;
			 
			}
	
		}
	
	}
			  
	return;
	 
}

void draw_line(double xi, double xf, double yi, double yf, int naxis1, int naxis2, int **data, int color, int thickness){

	/* Draws a line from (xi,yi) to (xf,yf), with a thickness of thickness. */
	
	int i,j;
	int x1,x2,y;
	double m,b;
	
	m=(yf-yi)/(xf-xi);
	b=yi-m*xi;
	
	x1=(int)xi;
	x2=(int)xf;
	
	for (i=x1;i<=x2;i++){
		
		/* Thicken the line a bit +/- thickness/2 pixels */
		
		y=(int)(m*(double)i+b);
		
		for(j=-1*(thickness/2);j<=(thickness/2);j++){
			
			if ((i > 0)&&(i < naxis1)&&((y+j) < naxis2)&&((y+j) > 0)){
			
				data[y+j][i]=color;
				
			}
		
		}	
		
	}
	
	return;	

}

void lint_to_binary(long number, int *binary){

	/* Converts a long integer to a binary number, returning an integer array
	 * of length 32 containing 0/1 depending upon whether or not a bit has 
	 * been set. */
	 
	int i;
	long quotient;
	 
	quotient=number;
	
	for (i=0;i<32;i++){
		
		binary[i]=0;
		 
		binary[i]=quotient%2;
		quotient=quotient/2;
		 
	}
	 
	return;
	
}

int good_pixel(long error_code){
	
	/* This function checks the array error_code (32 elements long) to see if a
	 * given pixel in SHARPs field data meets the following criteria: 1. It must
	 * be a good pixel (bit 27 is not set), 2. disambiguation must be either 
	 * trustworthy or doubtful (bit 6-7 either 00 or 01). 
	 * 
	 * See: http://hmi.stanford.edu/doc/magnetic/info_map_description.pdf	*/
	
	int *bin_error;
	
	bin_error=malloc(32*sizeof(int));
	
	lint_to_binary(error_code,bin_error);
	 
	if ((bin_error[27]==0)&&(bin_error[7]==0)){
		 
		/* 27:  not set, then a good pixel
		 * 7:	not set in trustworthy/doubtful cases */
		 
		free(bin_error);
		
		return 1;
		 
	} else {
		
		free(bin_error);
		
		return 0;
	
	}

}

double bilinear_interpolate(double x, double y, int naxis1, int naxis2, double **image){
	
	/* Simple bilinear interpolation routine. Keep in mind here that they are using a 
	 * (x,y) datapoint convention with the smallest being (1,1) and not (0,0), also 
	 * not [y][x] as used in image arrays. 	*/
	
	double f_R[2];
	double f_Q[2][2];
	double x1,x2,y1,y2;
	int xi,yi;
	double new_value;
	
	xi=(int)floor(x);
	yi=(int)floor(y);
	
	x1=(double)xi;
	x2=x1+1;
	y1=(double)yi;
	y2=y1+1;
	
	if ((xi == naxis1)||(xi == 0)||(yi == naxis2)||(yi == 0)){
		
		printf("Interpolation boundary issues.\n");
		
		new_value = image[yi][xi];

	} else {
		
		/* Using [x][y] for the f_Q's to keep with derivation */
		
		f_Q[0][0]=image[yi][xi];
		f_Q[1][0]=image[yi][xi+1];
		f_Q[0][1]=image[yi+1][xi];
		f_Q[1][1]=image[yi+1][xi+1];
		
		f_R[0]=((x2-x)/(x2-x1))*f_Q[0][0]+((x-x1)/(x2-x1))*f_Q[1][0];
		f_R[1]=((x2-x)/(x2-x1))*f_Q[0][1]+((x-x1)/(x2-x1))*f_Q[1][1];
		
		new_value=((y2-y)/(y2-y1))*f_R[0]+((y-y1)/(y2-y1))*f_R[1];
		
	}
	
	return new_value;
	
}

	


