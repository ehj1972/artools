/* Determines the area of the umbra and penumbra of both polarities within
 * a sunspot group.
 * 
 * 
 * Changelog:		
 * 			
 * 			2014.06.14		--	First version of this code.
 * 
 * 			2014.08.09		--		Added options to change sunspot cutoffs
 * 
 *
 * 			
 */

#define _GNU_SOURCE
#define image_bytes 1				/* One table of data as opposed to 3 for RGB */
#define colorscheme JCS_GRAYSCALE	/* One table = grayscale */
#define T_OBS_LENGTH 32				/* Character string length of T_OBS keyword */

#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <unistd.h>
#include <dirent.h>
#include <jpeglib.h>
#include "solarlib.h"

char *module_name="sharp_area";
char *module_desc="determines the area sunspots within an active region by using SHARPs data. Outputs data to ARnumber.area.dat.";
char *version_id="1.0";

ModuleArgs_t module_args[]={

    {ARG_STRING,"harpnum","","SHARPs number for the active region."},
	{ARG_STRING,"time_range","[]","time range to generate footpoints over (optional), format is [01.02.2013_04:05:06/??]"},
	{ARG_DOUBLE,"resolution","0.02058038057","latitude and longitude pixel resolution, in degrees. Default resolves to 250 km on a side. NOTE: can significantly increase run time."}, 
	{ARG_FLAG,"j","","When set, outputs projected images of continuum data. Image filename format is nnnnn.jpeg, where nnnnn is numerical ordering, starting with 00001.jpg. NOTE: higher resolution slows code considerably."},
	{ARG_DOUBLE,"umb_sig","2.575829","Standard error to base the cutoff between the quiet sun and what is considered to be a sunspot."},
	{ARG_DOUBLE,"pen_sig","1.0","Standard error used to determine what is considered to be an umbra and what is a penumbra."},
    {}

};

/* Main program */

int DoIt(){
 
    /* Variable declarations */

    /* Generic use variables */

    int i,j,k;					/* Generic loop counters */

    /* DRMS/JSOC variables */

	char *harpnum,*time_range;				/* SHARPs ID and time range from module_args[] */
	char *drms_query;							/* The DRMS query we are going to send */
	double dpix;								/* Resolution, in degrees, of a side of a pixel */
	double umb_sig, pen_sig;					/* Cutoffs for sunspots */
	int num_query_chars;						/* Number of characters in total query for malloc definition */
	CmdParams_t *params=&cmdparams;				/* For the module_args stuff */
	DRMS_RecordSet_t *drms_ids;					/* Holds the id structure of the records */
	DRMS_Record_t *drms_record;					/* Holds specific record information */
	DRMS_Segment_t *drms_segment;				/* Segment information */
	int drms_status;							/* Status variable */
	int num_records;							/* Number of records and variables to get the segment locations we need */
	int make_jpegs;								/* Flag for jpeg creation */
	
	/* cfitsio-related variables */

    fitsfile *c_ptr,*i_ptr;				/* Pointers for FITS files */
    int c_status, i_status;				/* Status checker for opening files */
    int hdu_pos,num_hdus;				/* HDU number and location */
    int needed_hdu;						/* The HDU we want to be on */
    int any_nulls;						/* Nonzero if there are null values in data */
    char i_filename[DRMS_MAXPATHLEN+1],c_filename[DRMS_MAXPATHLEN+1];	/* Containers for the filenames */
	double nulval;							/* Container for what null values in data should be set to */
	

    /* Needed FITS keywords */

    int naxis1,naxis2,naxis;			/* Length of axes */
    double crpix1,crpix2;				/* SHARPs CRPIX is the center of the Sun relative to the lower left of the patch (fortran: 1,1) */
    double imcrpix1,imcrpix2;			/* SHARPs IMCRPIX is the center of the Sun in full plate coordinates (fortran) */
    double cdelt1,cdelt2;				/* HMI pixel sizes in as */
    double crota2,crln_obs,crlt_obs;	/* Rotation, carrington observer location */
    double dsun_obs,rsun_obs;			/* Distance to and radius of the sun, observed at t_obs */
    double blank;						/* Values of bad data */
    TIME t_obs;							/* Observation time, ZULU */
    char **t_obs_s;						/* String version of observation time */
    long *quality;						/* QUALITY keyword */
    
    /* Data-related variables */
    
    double **con_dat, **inc_dat;			/* Data arrays, size defined once naxis1/2 are known */
    double *c_pixel_strip;					/* Strip of pixels read from the continuum data */
    double *i_pixel_strip;					/* Strip of pixels read from the inclination data */
    long fpixel[2];							/* Holds the location of pixel to be read */
    LONGLONG nelements;						/* Number of elements in a strip to be read (naxis2) */
    FILE *outptr,*bad_outptr;				/* Pointer to the outfile */
    char outfile_name[28];					/* Name of the outfile */
    char bad_outfile_name[24];				/* Name of the outfile for bad quality*/    
    
    /* Variables related to the actual task */
    
    double rllx,rlly;				/* Disk center, patch coordinates */
	double rscx,rscy;				/* Disk center, plate coordinates */
	double rpllx,rplly;				/* Lower left of patch, plate coordinates */
    double **radius;				/* Array to hold where the Sun's boundaries are, column major order */
    double rsun_pix,xx,yy,r;		/* Variables associated with radial boundaries */
    double tx,ty;					/* Intermediate coordinate conversion variables */
    double lat,lon;					/* Latitude and longitude variables */
    double min_lat,max_lat;			/* Min and max latitude, used to place boundaries on projection */
    double min_lon,max_lon;			/* Min and max longitude, used to place boundaries on projection */
    double num_good;				/* The number of nonzero continuum values in a sunspot region */
    double max_con_value;			/* Maximum continuum intensity value for jpeg creation */
    double umb,pen;					/* Intensity cutoffs for the umbra and penumbra */
    double *p_pen_area,*n_pen_area;	/* Area of polarity of penumbra within sunspot region */
    double *p_umb_area,*n_umb_area;	/* Area of polarity of umbra within sunspot region */
    double mu;						/* Angle from solar center */
    double x,y;						/* Placeholders for coordinates */
    int patch_x,patch_y;			/* New patch coordinates */
    int width,height;				/* Calculated array sizes for reprojection */
    long *t_secs;					/* T_OBS time converted to seconds since UNIX epoch */
    int noaa_ar;					/* Active region number */

    /* JPEG-related variables */
    
    FILE *jpeg_ptr;
    int **jpeg_outdata;
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    JSAMPROW row_pointer[1];
    char jpg_outfile[10];

    /* Initialise variables that need it */

    i=j=k=0;
    drms_status=0;

	/* Grab values from module_args[] */

	harpnum=strdup(params_get_str(params,"harpnum"));
	dpix=params_get_double(params,"resolution");
	time_range=strdup(params_get_str(params,"time_range"));
	make_jpegs=params_isflagset(params,"j");
	umb_sig=params_get_double(params,"umb_sig");
	pen_sig=params_get_double(params,"pen_sig");

	/* Now we start the main program. The rough layout is as follows: 
	   
	   1. Query drms. Based on what we learn, either exit or grab what we need from it.
	   2. Loop through the files, determining the number of sunspots in each polarity.
	   3. Print out a table of values at the end.				   */
	
	/* Forge the DRMS query. */
	
	num_query_chars=16+strlen(harpnum);		/* hmi.sharp_720s[] is 16 characters */
	
	if (time_range!="[]"){
		
		num_query_chars+=strlen(time_range);
		
	}
	
	drms_query=calloc(num_query_chars,sizeof(char));
	strcat(drms_query,"hmi.sharp_720s[");
	strcat(drms_query,harpnum);
	strcat(drms_query,"]");
	
	if (time_range!="[]"){
	
		strcat(drms_query,time_range);
		
	}
	
	/* Query DRMS for the records. */
	
	printf("Querying DRMS. This may take some time.\n");
	
	if (!(drms_ids=drms_open_records(drms_env,drms_query,&drms_status))){
		
		printf("No record sets match your criteria. Check your SHARPs ID and your time range.\n");
		drms_close_records(drms_ids,DRMS_FREE_RECORD);
		free(drms_query);
		return 0;
		
	}

	free(drms_query);

	num_records=drms_ids->n;
	printf("%d records match. Checking to make sure we have the needed data segments.\n",num_records);

	drms_record=drms_ids->records[0];

	if (!(drms_segment_lookup(drms_record,"continuum"))){
		
		printf("Continuum segment not present! Exiting.\n");
		drms_close_records(drms_ids,DRMS_FREE_RECORD);
		return 0;
	
	}

	if (!(drms_segment_lookup(drms_record,"inclination"))){
		
		printf("Inclination segment not present! Exiting.\n");
		drms_close_records(drms_ids,DRMS_FREE_RECORD);
		return 0;
	
	}
	
	/* Now we can start. We loop over all records, get the keyword information we need, get the data, and then calculate. */
	
	/* First set the size of some of the pointers that are intended to be arrays. */
	
	/* Data-related arrays */

	p_pen_area=malloc(num_records*sizeof(double));
	n_pen_area=malloc(num_records*sizeof(double));
	
	p_umb_area=malloc(num_records*sizeof(double));
	n_umb_area=malloc(num_records*sizeof(double));	
	
	/* Time-related arrays */
	
	t_secs=malloc(num_records*sizeof(long));
	t_obs_s=malloc(num_records*sizeof(char *));
	
	for (i=0;i<num_records;i++){
		
		t_obs_s[i]=malloc(T_OBS_LENGTH*sizeof(char));
		
	}
	
	/* Quality */
	
	quality=malloc(num_records*sizeof(long));
		
	for (k=0;k<num_records;k++){ 

		drms_record=drms_ids->records[k];
		
		/* Check to see if we have data for a particular observation time */
		
		if (!(t_obs=drms_getkey_time(drms_record,"T_OBS",&drms_status))){
					
			printf("Keyword %s not present! Exiting.\n","T_OBS");
			drms_close_records(drms_ids,DRMS_FREE_RECORD);
			return 0;
		
		}
		
		if (time_is_invalid(t_obs)){
			
			printf("Bad record, skipping.\n");
			
		} else {

			if (k == 0){
				
				/* Get active region number for pregame info */
			
				if (!(noaa_ar=drms_getkey_int(drms_record,"NOAA_AR",&drms_status))){
					
					printf("Keyword %s not present! Exiting.\n","NOAA_AR");
					drms_close_records(drms_ids,DRMS_FREE_RECORD);
					return 0;
		
				}
				
				/* Create data output filename */
			
				sprintf(outfile_name,"AR%d.area-output.dat",noaa_ar);
				sprintf(bad_outfile_name,"AR%d.bad-quality.dat",noaa_ar);
				
				/* Print out some stats before we run */
			
				printf("Pre-analysis details:\n");
				printf("Analyzing: HARP Number: %s\t NOAA AR: %d\n",harpnum,noaa_ar);
				printf("Output file: %s\n",outfile_name);
			
				if (time_range == "[]"){
				
					printf("We are covering the entire data set.\n");
			
				} else {
				
					printf("Time range: %s\n",time_range);
			
				}

					printf("\n\n");
			
			}
			
			printf("Processing record %d of %d.\n",k+1,num_records);
		
			/* Fill keywords */
			
			if (!(quality[k]=drms_getkey_longlong(drms_record,"QUALITY",&drms_status))){
			
				quality[k]=0;
		
			}
	
			if (!(crpix1=drms_getkey_double(drms_record,"CRPIX1",&drms_status))){
			
				printf("Keyword %s not present! Exiting.\n","CRPIX1");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
			
			if (!(crpix2=drms_getkey_double(drms_record,"CRPIX2",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","CRPIX2");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}

			if (!(imcrpix1=drms_getkey_double(drms_record,"IMCRPIX1",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","IMCRPIX1");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
		
			if (!(imcrpix2=drms_getkey_double(drms_record,"IMCRPIX2",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","IMCRPIX2");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}

		
			if (!(cdelt1=drms_getkey_double(drms_record,"CDELT1",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","CDELT1");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
			
			if (!(cdelt2=drms_getkey_double(drms_record,"CDELT2",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","CDELT2");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}

			if (!(crota2=drms_getkey_double(drms_record,"CROTA2",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","CROTA2");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
		
			if (!(crln_obs=drms_getkey_double(drms_record,"CRLN_OBS",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","CRLN_OBS");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
		
			if (!(crlt_obs=drms_getkey_double(drms_record,"CRLT_OBS",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","CRLT_OBS");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
		
			if (!(rsun_obs=drms_getkey_double(drms_record,"RSUN_OBS",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","RSUN_OBS");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
			
			if (!(dsun_obs=drms_getkey_double(drms_record,"DSUN_OBS",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","DSUN_OBS");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
		
			if (!(blank=drms_getkey_double(drms_record,"BLANK",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","BLANK");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
		
			sprint_ut(t_obs_s[k],t_obs);		/* Convert to string */
		
			printf(" Checking quality.\n");
			
			if (quality[k] == 0){
		
				/* Ok, now we get the file locations for the data from DRMS. Originally I attempted
				* to pull the data directly from DRMS but there were vague free() and malloc() issues
				* possibly due to icc and multithreading. */
		
				if (!(drms_segment=drms_segment_lookup(drms_record,"continuum"))){
			
					printf("Problem opening the continuum segment! Exiting.\n");
					drms_close_records(drms_ids,DRMS_FREE_RECORD);
					return 0;
			
				}
		
				drms_segment_filename(drms_segment,c_filename);
			
				/* naxis1 and naxis2 we get from segment information */
		
				naxis1=(int)drms_segment->axis[0];
				naxis2=(int)drms_segment->axis[1];
		
				if (!(drms_segment=drms_segment_lookup(drms_record,"inclination"))){
			
					printf("Problem opening the inclination segment! Exiting.\n");
					drms_close_records(drms_ids,DRMS_FREE_RECORD);
					return 0;
			
				}
		
				drms_segment_filename(drms_segment,i_filename);
			
				printf(" Opening FITs files.\n");
		
				/* Now open the FITS files and get to work. */
	
				c_status=0;
				i_status=0;

				if (fits_open_file(&c_ptr,c_filename,READONLY,&c_status)){
				
					printf("Cannot open %s! Exiting.\n",c_filename);
					exit(0);

				}

				if (fits_open_file(&i_ptr,i_filename,READONLY,&i_status)){
				
					printf("Cannot open %s! Exiting.\n",i_filename);
					fits_close_file(c_ptr,&c_status);
					exit(0);

				}


				/* Walk through the HDUs until we get to one with NAXIS=2 */

				/* Number of headers and which one we are on */

				fits_get_num_hdus(c_ptr,&num_hdus,&c_status);
				fits_get_hdu_num(c_ptr,&hdu_pos);

				/* Find the one with two naxes */

				naxis=0;
			
				for (i=1;i<=num_hdus;i++){
				
					fits_movabs_hdu(c_ptr,i,NULL,&c_status);
					fits_read_key(c_ptr,TINT,"NAXIS",&naxis,NULL,&c_status);
				
					if (naxis == 2){
	
						needed_hdu=i;
	
					}
				}

				/* Set all to the needed HDU */

				if (naxis == 0){
		
					printf("HDU problems: can't find one with the required number of axes.\n");
					fits_close_file(c_ptr,&c_status);
					fits_close_file(i_ptr,&i_status);
					exit(0);

				} else {

					fits_movabs_hdu(c_ptr,needed_hdu,NULL,&c_status);
					fits_movabs_hdu(i_ptr,needed_hdu,NULL,&i_status);

				}

				/* Now set some definite boundaries on arrays to help with memory usage. Arrays are called as con_dat[y][x]. */
	
				con_dat=malloc(naxis2*sizeof(double *));
				inc_dat=malloc(naxis2*sizeof(double *));
				radius=malloc(naxis2*sizeof(double *));
				
				c_pixel_strip=calloc(naxis1,sizeof(double));
				i_pixel_strip=calloc(naxis1,sizeof(double));
	
				for (j=0;j<naxis2;j++){
		
					con_dat[j]=calloc(naxis1,sizeof(double));
					inc_dat[j]=calloc(naxis1,sizeof(double));
					radius[j]=calloc(naxis1,sizeof(double));
		
				}
    
				/* Next, put data in the arrays, clear out blank values and NaNs from the data */
	
				fpixel[0]=1;		/* cfitsio uses fortran reference instead of c when accessing data */
				nelements=naxis1;
				nulval=0;
				any_nulls=0;
    
				for (j=0;j<naxis2;j++){
		
					fpixel[1]=j+1;	/* Add 1 to account for fortran FITS indexing */
					fits_read_pix(c_ptr,TDOUBLE,fpixel,nelements,&nulval,c_pixel_strip,&any_nulls,&c_status);
					fits_read_pix(i_ptr,TDOUBLE,fpixel,nelements,&nulval,i_pixel_strip,&any_nulls,&i_status);
		
					for (i=0;i<naxis1;i++){
			
						/* Kill blank values if present, if not assign to the correct place in array */
			
						if (c_pixel_strip[i] == blank){
				
							con_dat[j][i]=0.0;
					
						} else {
					
							con_dat[j][i]=c_pixel_strip[i];
			
						}
			
						if (i_pixel_strip[i] == blank){
				
							inc_dat[j][i]=0.0;
				
						} else {
					
							inc_dat[j][i]=i_pixel_strip[i];
			
						}
					
					}
				}
		
				printf(" Analyzing data.\n");
	
				/* We now have filled arrays of data about the SHARPs patch and all of the needed keywords. */
	
				/* Coordinate variable definitions for future use. See the declaration section */
				/* at the front of the code to understand what these refer to. */
	
				rscx=(imcrpix1-1);
				rscy=(imcrpix2-1);
				rllx=(crpix1-1);
				rlly=(crpix2-1);
				rpllx=(rscx-rllx);
				rplly=(rscy-rlly);
	
				/* Fill the radius array, handle limb darkening, find max value in patch. */
	
				rsun_pix=rsun_obs/cdelt1;
				max_con_value=0.0;
	
				for (i=0;i<naxis1;i++){
		
					for (j=0;j<naxis2;j++){
			
						xx=(double)i-rllx;
						yy=(double)j-rlly;
						r=sqrt(xx*xx+yy*yy);
				
						if (r <= rsun_pix){
				
							radius[j][i]=r;
							mu=sqrt(1-pow(r/rsun_pix,2.0));
							con_dat[j][i]=con_dat[j][i]/limb_darken(mu);

							if (con_dat[j][i] > max_con_value){
				
								max_con_value=con_dat[j][i];
				
							}
					
						}
					}
				}
				
				/* Next, set latitude and longitude boundaries on the patch */
				
				min_lat=1000;
				max_lat=-1000;
				min_lon=1000;
				max_lon=-1000;
	
				for (i=0;i<naxis1;i++){
		
					for (j=0;j<naxis2;j++){
			
						x=rpllx+(double)i;
						y=rplly+(double)j;
			
						pixel2hpl(x,y,crota2,imcrpix1,imcrpix2,cdelt1,cdelt2,&tx,&ty);
						hpl2stony(tx,ty,dsun_obs,crln_obs,crlt_obs,&lat,&lon);
						
						lon-=crln_obs;
						
						/* Lat and lon -1000's come from nonphysical solutions */
			
						if ((lat < min_lat)&&(lat > -360)){
				
							min_lat=lat;
			
						}
			
						if ((lon < min_lon)&&(lon > -360)){
				
							min_lon=lon;
			
						}
			
						if ((lat > max_lat)&&(lat > -360)){
				
							max_lat=lat;
				
						}
			
						if ((lon > max_lon)&&(lon > -360)){
				
							max_lon=lon;
				
						}
			
					}
				}
				
				/* New width and heights for area counters and -j flag, padding by a pixel on each side */
				
				width=(int)rint((max_lon-min_lon)/dpix)+2;
				height=(int)rint((max_lat-min_lat)/dpix)+2;
				
				/* Create image array if -j is set */
				
				if (make_jpegs){
					
					row_pointer[0]=(unsigned char*)malloc(width*image_bytes);
					
					jpeg_outdata=malloc(height*sizeof(int *));
					
					for (i=0;i<height;i++){
						
						jpeg_outdata[i]=calloc(width,sizeof(int));
						
					}
				
				}
					
				/* Create the actual masks for sunspot counting */
			
				spotbounds(con_dat,naxis1,naxis2,&umb,&pen,umb_sig,pen_sig);

				num_good=0;
	
				for (i=0;i<naxis1;i++){
		
					for (j=0;j<naxis2;j++){
		
						if ((radius[j][i] > 0)&&(con_dat[j][i] <= pen)&&(con_dat[j][i] > 0)){
				
							num_good++;
			
						}
					}
				}
				
				p_pen_area[k]=0;
				n_pen_area[k]=0;
				p_umb_area[k]=0;
				n_umb_area[k]=0;
	
				if (num_good > 0){
		
					/* We have a region. */
					
					/* Now, the area. We are not looping over the image per se, we are looping
					 * over the projection, looking up a lat and lon, and then seeing if that lat and lon
					 * lives in the image. We round past 0.5 pixels to the next pixel, so there isn't 
					 * any interpolation going on. The check as to whether or not it is counted is 
					 * whether or not the pixel, once it passes the aformentioned check, also lives
					 * within the umbra or the penumbra. Area is calculated as a spherical segment, 
					 * e.g. dA=R^2\cos(theta)d\thetad\phi */
					 
					 
					for (i=1;i<(width-1);i++){
		
						for (j=1;j<(height-1);j++){
			
							/* Ok, x runs from min_lon to max_lon, and y runs from min_lat to max_lat. Let's 
							* cut a touch from the boundaries so we don't get wierd boundary issues. */
			 
							lon=((double)i*dpix)+min_lon+crln_obs;	
							lat=((double)j*dpix)+min_lat;	
			
							stony2hpl(lat,lon,dsun_obs,crln_obs,crlt_obs,&tx,&ty);
							hpl2pixel(tx,ty,crota2,imcrpix1,imcrpix2,cdelt1,cdelt2,&x,&y);
			
							patch_x=(int)rint(x-rpllx);
							patch_y=(int)rint(y-rplly);
							
							/* Zero jpeg stuff */
							
							if (make_jpegs){
								
								jpeg_outdata[j][i]=0;
								
							}
			
							/* Got the patch x,y locations. Now sanity check to see if they actually live within the patch
							 * and are on the Sun. */ 
							
							if ((patch_x < naxis1)&&(patch_x >= 0)&&(patch_y < naxis2)&&(patch_y >= 0)){
								
								/* Now, find whether the pixel in question is umbra, penumbra, or background */
								
								if (make_jpegs){
				
									jpeg_outdata[j][i]=(int)(256*con_dat[patch_y][patch_x]/max_con_value);
									
								}
								
								if (radius[patch_y][patch_x] < rsun_pix){
								
									if ((con_dat[patch_y][patch_x] <= pen)&&(con_dat[patch_y][patch_x] > umb)){
									
										if (inc_dat[patch_y][patch_x] < 90){
									
											/* "Positive" penumbra pixel */
									
											p_pen_area[k]+=pixel_area(lat,lon,dpix);
									
										}
									
										if (inc_dat[patch_y][patch_x] > 90){
									
											/* "Negative" penumbra pixel */
									
											n_pen_area[k]+=pixel_area(lat,lon,dpix);
									
										}
								
									}								
	
									if ((con_dat[patch_y][patch_x] <= umb)&&(con_dat[patch_y][patch_x] > 0)){
									
										if (inc_dat[patch_y][patch_x] < 90){
									
											/* "Positive" umbra pixel */
									
											p_umb_area[k]+=pixel_area(lat,lon,dpix);
									
										}
									
										if (inc_dat[patch_y][patch_x] > 90){
									
											/* "Negative" umbra pixel */
									
											n_umb_area[k]+=pixel_area(lat,lon,dpix);
									
										}
								
									}
									
								}	
				
							} 
			
						}
					}
					 
				}

				if (make_jpegs){
					
					printf(" Creating jpeg image\n");
				
					sprintf(jpg_outfile,"%05d.jpg",k+1);
			
					jpeg_ptr=fopen(jpg_outfile,"wb");
	
					cinfo.err=jpeg_std_error(&jerr);
					jpeg_create_compress(&cinfo);
					jpeg_stdio_dest(&cinfo,jpeg_ptr);
	
					cinfo.image_width=width;
					cinfo.image_height=height;
					cinfo.input_components=image_bytes;
					cinfo.in_color_space=colorscheme;
	
					jpeg_set_defaults(&cinfo);
					jpeg_set_quality(&cinfo,80,TRUE);
					jpeg_start_compress(&cinfo,TRUE);
	
					for (j=0;j<height;j++){
						
						for (i=0;i<width;i++){
							
							row_pointer[0][i]=(unsigned char)jpeg_outdata[j][i];
						
						}
		
						jpeg_write_scanlines(&cinfo,row_pointer,1);
		
					}
	
					jpeg_finish_compress(&cinfo);
					jpeg_destroy_compress(&cinfo);
				
					fclose(jpeg_ptr);
				}
			
				/* Free up dynamic arrays and close the files for the next iteration. */
				
				for (i=0;i<naxis2;i++){
				
					free(con_dat[i]);
					free(inc_dat[i]);
					free(radius[i]);
					
					if (make_jpegs){
					
						free(jpeg_outdata[i]);
				
					}
					
					
					
				
				}
				
				if (make_jpegs){
					
					free(jpeg_outdata);
					free(row_pointer[0]);
				
				}
				
				free(con_dat);
				free(inc_dat);
				free(radius);
				
				free(c_pixel_strip);
				free(i_pixel_strip);

				fits_close_file(c_ptr,&c_status);
				fits_close_file(i_ptr,&i_status);

			}
	
		}
		printf(" Completed run %d of %d.\n",k+1,num_records);
	}
	
	/* Output the results to a text file. */
	
	/* First, convert t_obs to a more reasonable product. */
	
	for (i=0;i<num_records;i++){			
	
		t_secs[i]=isotime(t_obs_s[i]);

	}

	printf("Writing to file %s.\n",outfile_name);
	
	outptr=fopen(outfile_name,"w");
	bad_outptr=fopen(bad_outfile_name,"w");
	
	fprintf(outptr,"Active Region: %d NOTE: Areas in km^2\n",noaa_ar);
	fprintf(outptr,"Full UT Time            time(s)  p umb area  n umb area  p pen area  n pen area\n");
	fprintf(outptr,"----------------------  -------  ----------  ----------  ----------  ----------\n");
	
	fprintf(bad_outptr,"Active Region: AR%d\n",noaa_ar);
	fprintf(bad_outptr,"Full UT Time            quality  time(s)  number\n");
	fprintf(bad_outptr,"----------------------  -------  -------  ------\n");
	
	for (i=0;i<num_records;i++){
		
		/* Make times relative to the first recorded one. */
		
		if (quality[i] == 0){
		
			fprintf(outptr,"%s  %7ld  %7.2lf  %7.2lf  %7.2lf  %7.2lf\n",t_obs_s[i],t_secs[i]-t_secs[0],p_umb_area[i]/1000000.0,n_umb_area[i]/1000000.0,p_pen_area[i]/1000000.0,n_pen_area[i]/1000000.0);
	
		} else {
			
			fprintf(bad_outptr,"%s  0x%08lX  %7ld  %d\n",t_obs_s[i],quality[i],t_secs[i],i+1);
			
		}
	}
	
	fclose(outptr);
	fclose(bad_outptr);
	
	printf("Done!\n");
	
	/* Free memory associated with dynamic arrays. */
    
	free(p_umb_area);
	free(n_umb_area);
	free(p_pen_area);
	free(n_pen_area);

	free(t_secs);
	free(t_obs_s);
	
	free(quality);
	
	/* Close the drms records connection */
	
	drms_close_records(drms_ids,DRMS_FREE_RECORD);
    
	/* Done! */
    
	return 0;
}
