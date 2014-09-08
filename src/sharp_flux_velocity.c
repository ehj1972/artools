/* This code determines what the average LOS velocities are of plasma within the umbral and penumbral regions of a sunspot
 * group, separated by preceding and follower. 
 * 
 * All code was originally written in Python/Numpy or IDL and the ported to C, and then tested 
 * accordingly. Compile and run sharps_total_field to see the command line options.
 * 
 * 
 * Changelog:		
 * 			
 * 			2014.09.03		--		First version of this code.	
 * 		
 */

#define _GNU_SOURCE
#define T_OBS_LENGTH 32				/* Character string length of T_OBS keyword */

#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <unistd.h>
#include <dirent.h>
#include "solarlib.h"

char *module_name="sharp_field";
char *module_desc="Calculates various field products of both polarities within a sunspot group. Outputs data to ARnumber.field.dat.";
char *version_id="1.0";

ModuleArgs_t module_args[]={

    {ARG_STRING,"harpnum","","SHARPs number for the active region."},
	{ARG_STRING,"time_range","[]","time range to generate footpoints over (optional), format is [01.02.2013_04:05:06/??]"},
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

	char *harpnum,*time_range;					/* SHARPs ID and time range from module_args[] */
	double umb_sig, pen_sig;					/* Cutoffs for sunspots */
	char *drms_query;							/* The DRMS query we are going to send */
	int num_query_chars;						/* Number of characters in total query for malloc definition */
	CmdParams_t *params=&cmdparams;				/* For the module_args stuff */
	DRMS_RecordSet_t *drms_ids;					/* Holds the id structure of the records */
	DRMS_Record_t *drms_record;					/* Holds specific record information */
	DRMS_Segment_t *drms_segment;				/* Segment information */
	int drms_status;							/* Status variable */
	int num_records;							/* Number of records and variables to get the segment locations we need */
	
	/* cfitsio-related variables */

    fitsfile *v_ptr,*i_ptr,*c_ptr;		/* Pointers for FITS files */
    fitsfile *info_ptr;		
    int v_status,i_status,c_status;		/* Status checker for opening files */
    int info_status;					
    int hdu_pos,num_hdus;				/* HDU number and location */
    int needed_hdu;						/* The HDU we want to be on */
    int any_nulls;						/* Nonzero if there are null values in data */
    char i_filename[DRMS_MAXPATHLEN+1],v_filename[DRMS_MAXPATHLEN+1];		/* Containers for the filenames */
	char c_filename[DRMS_MAXPATHLEN+1],info_filename[DRMS_MAXPATHLEN+1];
	double nulval;							/* Container for what null values in data should be set to */
	long l_nulval;

    /* Needed FITS keywords */

    int naxis1,naxis2,naxis;			/* Length of axes */
    double crpix1,crpix2;				/* SHARPs CRPIX is the center of the Sun relative to the lower left of the patch (fortranL 1,1) */
    double imcrpix1,imcrpix2;			/* SHARPs IMCRPIX is the center of the Sun in full plate coordinates (fortran) */
    double cdelt1,cdelt2;				/* HMI pixel sizes in as */    double blank;						/* Values of bad data */
    double crota2,crln_obs,crlt_obs;	/* Rotation, carrington observer location */
    double dsun_obs,rsun_obs;			/* Distance to and radius of the sun, observed at t_obs */
    double obs_vr,obs_vw,obs_vn;		/* Satellite velocity keywords */
    TIME t_obs;							/* Observation time, ZULU */
    char **t_obs_s;						/* String version of observation time */
    long *quality;						/* Data quality keyword */
    
    /* Data-related variables */
    
    double **vlos_dat, **inc_dat;			/* Data arrays, size defined once naxis1/2 are known */
    double **con_dat;
    long **info_dat;
    double *v_pixel_strip;					/* Strip of pixels read from the field data */
    double *i_pixel_strip;					/* Strip of pixels read from the inclination data */
    double *c_pixel_strip;					/* Strip of pixels read from the continuum data */
    long *info_pixel_strip;					/* Strip of pixels read from the continuum data */
    long fpixel[2];							/* Holds the location of pixel to be read */
    LONGLONG nelements;						/* Number of elements in a strip to be read (naxis2) */
    FILE *outptr,*bad_outptr;				/* Pointer to the outfile */
    char outfile_name[25];					/* Name of the outfile */
    char bad_outfile_name[24];				/* Name of the outfile for bad quality*/
    
    /* Variables related to the actual task */
    
    double rscx,rscy;				/* Disk center, plate coordinates */
    double rllx,rlly;				/* Disk center, patch coordinates */
    double rpllx,rplly;				/* Lower left of patch, plate coordinates */
    double **radius;				/* Array to hold where the Sun's boundaries are, column major order */
	double umb,pen;					/* Intensity cutoffs for the umbra and penumbra */
    double num_good;				/* The number of nonzero continuum values in a sunspot region */
    double mu;						/* Angle from solar center */
    double rsun_pix,xx,yy,r;		/* Variables associated with radial boundaries */
    long *t_secs;					/* T_OBS time converted to seconds since UNIX epoch */
    double *umb_p_vel,*umb_n_vel;	/* Umbral LOS velocities */
    double *pen_p_vel,*pen_n_vel;	/* Penumbral LOS velocities */
    double *f_umb_p_vel,*f_umb_n_vel;	/* Umbral LOS velocities, "fixed" */
    double *f_pen_p_vel,*f_pen_n_vel;	/* Penumbral LOS velocities, "fixed" */
    double umb_p_counter,umb_n_counter;	/* Counters for averages */
    double pen_p_counter,pen_n_counter;
    double vr,vn,vw,vt;				/* Various LOS velocities due to satellite motion and solar rotation */
    double vlos;					/* VLOS after removing all of the satellite motion and solar rotation */
    double vlos_fixed;				/* VLOS after removing all of the satellite motion and solar rotation and flattened using mu position */
    double x,y;						/* Pixel location on the solar disk */
    double ttx,tty,tx,ty,zz;		/* Intermediate variables for coordinate conversions */
    double lat,lon;					/* Carrington coordinate latitude and longitude */
    double p,t;						/* Spherical phi and theta coordinates */
	int noaa_ar;					/* Active region number */
       
    /* Initialise variables that need it */

    i=j=k=0;
    drms_status=0;
    
	/* Grab values from module_args[] */

	harpnum=strdup(params_get_str(params,"harpnum"));
	time_range=strdup(params_get_str(params,"time_range"));
	umb_sig=params_get_double(params,"umb_sig");
	pen_sig=params_get_double(params,"pen_sig");
	
	/* Now we start the main program. The rough layout is as follows: 
	   
	   1. Query drms. Based on what we learn, either exit or grab what we need from it.
	   2. Loop through the files, calculate field stuff.
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

	if (!(drms_segment_lookup(drms_record,"vlos_mag"))){
		
		printf("vlos_mag segment not present! Exiting.\n");
		drms_close_records(drms_ids,DRMS_FREE_RECORD);
		return 0;
	
	}

	if (!(drms_segment_lookup(drms_record,"inclination"))){
		
		printf("Inclination segment not present! Exiting.\n");
		drms_close_records(drms_ids,DRMS_FREE_RECORD);
		return 0;
	
	}

	if (!(drms_segment_lookup(drms_record,"continuum"))){
		
		printf("Continuum segment not present! Exiting.\n");
		drms_close_records(drms_ids,DRMS_FREE_RECORD);
		return 0;
	
	}

	if (!(drms_segment_lookup(drms_record,"info_map"))){
		
		printf("Info_map segment not present! Exiting.\n");
		drms_close_records(drms_ids,DRMS_FREE_RECORD);
		return 0;

	}

	/* Now we can start. We loop over all records, get the keyword information we need, get the data, and then calculate. */
	
	/* First set the size of some of the pointers that are intended to be arrays. */
	
	/* Time-related arrays */
	
	t_secs=malloc(num_records*sizeof(long));
	t_obs_s=malloc(num_records*sizeof(char *));
	
	for (i=0;i<num_records;i++){
		
		t_obs_s[i]=malloc(T_OBS_LENGTH*sizeof(char));
		
	}
	
	/* Data-related arrays */
	
	umb_p_vel=malloc(num_records*sizeof(double));
	pen_p_vel=malloc(num_records*sizeof(double));
	umb_n_vel=malloc(num_records*sizeof(double));
	pen_n_vel=malloc(num_records*sizeof(double));
	
	f_umb_p_vel=malloc(num_records*sizeof(double));
	f_pen_p_vel=malloc(num_records*sizeof(double));
	f_umb_n_vel=malloc(num_records*sizeof(double));
	f_pen_n_vel=malloc(num_records*sizeof(double));

	/* Record quality */
	
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
			
				sprintf(outfile_name,"AR%d.fvlos-output.dat",noaa_ar);
				sprintf(bad_outfile_name,"AR%d.bad-quality.dat",noaa_ar);
							
				/* Print out some stats before we run */
			
				printf("Pre-analysis details:\n");
				printf("Analyzing: HARP Number: %s\t NOAA AR: %d\n",harpnum,noaa_ar);
				printf("Output file: %s\n",outfile_name);
			
				if (time_range == "[]\0"){
				
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
			
			if (!(obs_vr=drms_getkey_double(drms_record,"OBS_VR",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","OBS_VR");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}

			if (!(obs_vn=drms_getkey_double(drms_record,"OBS_VN",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","OBS_VN");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
			
			if (!(obs_vw=drms_getkey_double(drms_record,"OBS_VW",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","OBS_VW");
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
				* possibly due to icc and multithreading, or library problems. */
		
				if (!(drms_segment=drms_segment_lookup(drms_record,"vlos_mag"))){
			
					printf("Problem opening the vlos_mag segment! Exiting.\n");
					drms_close_records(drms_ids,DRMS_FREE_RECORD);
					return 0;
			
				}
		
				drms_segment_filename(drms_segment,v_filename);
			
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
			
				if (!(drms_segment=drms_segment_lookup(drms_record,"info_map"))){
			
					printf("Problem opening the info_map segment! Exiting.\n");
					drms_close_records(drms_ids,DRMS_FREE_RECORD);
					return 0;
			
				}
			
				drms_segment_filename(drms_segment,info_filename);
			
				printf(" Opening FITs files.\n");
		
				/* Now open the FITS files and get to work. */
	
				v_status=0;
				i_status=0;
				c_status=0;
				info_status=0;
				any_nulls=0;

				if (fits_open_file(&v_ptr,v_filename,READONLY,&v_status)){
				
					printf("Cannot open %s! Exiting.\n",v_filename);
					exit(0);

				}

				if (fits_open_file(&i_ptr,i_filename,READONLY,&i_status)){
				
					printf("Cannot open %s! Exiting.\n",i_filename);
					fits_close_file(i_ptr,&i_status);
					exit(0);

				}
			
				if (fits_open_file(&c_ptr,c_filename,READONLY,&c_status)){
				
					printf("Cannot open %s! Exiting.\n",c_filename);
					fits_close_file(c_ptr,&c_status);
					exit(0);

				}
			
				if (fits_open_file(&info_ptr,info_filename,READONLY,&info_status)){
				
					printf("Cannot open %s! Exiting.\n",info_filename);
					fits_close_file(info_ptr,&info_status);
					exit(0);

				}
				
				
				/* Walk through the HDUs until we get to one with NAXIS=2 */

				/* Number of headers and which one we are on */

				fits_get_num_hdus(v_ptr,&num_hdus,&v_status);
				fits_get_hdu_num(v_ptr,&hdu_pos);

				/* Find the one with two naxes */

				naxis=0;

				for (i=1;i<=num_hdus;i++){
				
					fits_movabs_hdu(v_ptr,i,NULL,&v_status);
					fits_read_key(v_ptr,TINT,"NAXIS",&naxis,NULL,&v_status);
				
					if (naxis == 2){
	
						needed_hdu=i;
	
					}
				}


				/* Set all to the needed HDU */

				if (naxis == 0){
		
					printf("HDU problems: can't find one with the required number of axes.\n");
					fits_close_file(v_ptr,&v_status);
					fits_close_file(i_ptr,&i_status);
					fits_close_file(c_ptr,&c_status);
					fits_close_file(info_ptr,&info_status);
					exit(0);

				} else {

					fits_movabs_hdu(v_ptr,needed_hdu,NULL,&v_status);
					fits_movabs_hdu(i_ptr,needed_hdu,NULL,&i_status);
					fits_movabs_hdu(c_ptr,needed_hdu,NULL,&c_status);
					fits_movabs_hdu(info_ptr,needed_hdu,NULL,&info_status);

				}


				/* Now set some definite boundaries on arrays to help with memory usage. Arrays are called as con_dat[y][x]. */
				
				
				vlos_dat=malloc(naxis2*sizeof(double *));
				inc_dat=malloc(naxis2*sizeof(double *));
				con_dat=malloc(naxis2*sizeof(double *));
				radius=malloc(naxis2*sizeof(double *));
				info_dat=malloc(naxis2*sizeof(long *));
			
				v_pixel_strip=calloc(naxis1,sizeof(double));
				i_pixel_strip=calloc(naxis1,sizeof(double));
				c_pixel_strip=calloc(naxis1,sizeof(double));
				info_pixel_strip=calloc(naxis1,sizeof(long));
	
				for (j=0;j<naxis2;j++){
		
					vlos_dat[j]=calloc(naxis1,sizeof(double));
					inc_dat[j]=calloc(naxis1,sizeof(double));
					con_dat[j]=calloc(naxis1,sizeof(double));
					radius[j]=calloc(naxis1,sizeof(double));
					info_dat[j]=calloc(naxis1,sizeof(long));
		
				}
				
				
				/* Next, put data in the arrays, clear out blank values and NaNs from the data */
	
				fpixel[0]=1;		/* cfitsio uses fortran reference instead of c when accessing data */
				nelements=naxis1;
				nulval=0;
				l_nulval=0;
				any_nulls=0;
    
				for (j=0;j<naxis2;j++){
		
					fpixel[1]=j+1;	/* Add 1 to account for fortran FITS indexing */
					fits_read_pix(v_ptr,TDOUBLE,fpixel,nelements,&nulval,v_pixel_strip,&any_nulls,&v_status);
					fits_read_pix(i_ptr,TDOUBLE,fpixel,nelements,&nulval,i_pixel_strip,&any_nulls,&i_status);
					fits_read_pix(c_ptr,TDOUBLE,fpixel,nelements,&nulval,c_pixel_strip,&any_nulls,&c_status);
					fits_read_pix(info_ptr,TLONG,fpixel,nelements,&l_nulval,info_pixel_strip,&any_nulls,&info_status);
		
					for (i=0;i<naxis1;i++){
			
						/* Kill blank values if present, if not assign to the correct place in array */
			
						if (v_pixel_strip[i] == blank){
				
							vlos_dat[j][i]=0.0;
					
						} else {
				
							vlos_dat[j][i]=v_pixel_strip[i];
			
						}
			
						if (i_pixel_strip[i] == blank){
				
							inc_dat[j][i]=0.0;
				
						} else {
					
							inc_dat[j][i]=i_pixel_strip[i];
			
						}

						if (c_pixel_strip[i] == blank){
				
							con_dat[j][i]=0.0;
			
						} else {
					
							con_dat[j][i]=c_pixel_strip[i];
			
						}
					
						info_dat[j][i]=info_pixel_strip[i];
					
					}
				}
				
		
				printf(" Analyzing data.\n");
	
				/* Coordinate variable definitions for future use. See the declaration section */
				/* at the front of the code to understand what these refer to. */
	
				rscx=(imcrpix1-1.0);
				rscy=(imcrpix2-1.0);
				rllx=(crpix1-1.0);
				rlly=(crpix2-1.0);
				rpllx=(rscx-rllx);
				rplly=(rscy-rlly);
	
				/* Fill the radius array, handle limb darkening */
	
				rsun_pix=rsun_obs/cdelt1;
	
				for (i=0;i<naxis1;i++){
		
					for (j=0;j<naxis2;j++){
			
						xx=(double)i-rllx;
						yy=(double)j-rlly;
						r=sqrt(xx*xx+yy*yy);
				
						if (r <= rsun_pix){
				
							radius[j][i]=r;
							mu=sqrt(1-pow(r/rsun_pix,2.0));
							con_dat[j][i]=con_dat[j][i]/limb_darken(mu);
					
						}
					}
				}
			
				/* Determine if we actually have an active region here. */
			
				spotbounds(con_dat,naxis1,naxis2,&umb,&pen,umb_sig,pen_sig);
			
				num_good=0;
	
				for (i=0;i<naxis1;i++){
		
					for (j=0;j<naxis2;j++){
		
						if ((radius[j][i] > 0)&&(con_dat[j][i] <= pen)&&(con_dat[j][i] > 0)){
				
							num_good++;
			
						}
					}
				}
	
				/* Now find the average LOS velocities. */
			
				umb_p_vel[k]=0;
				pen_p_vel[k]=0;
				umb_n_vel[k]=0;
				pen_n_vel[k]=0;
				
				f_umb_p_vel[k]=0;
				f_pen_p_vel[k]=0;
				f_umb_n_vel[k]=0;
				f_pen_n_vel[k]=0;
				
				umb_p_counter=0.0;
				umb_n_counter=0.0;
				pen_p_counter=0.0;
				pen_n_counter=0.0;
			
				if (num_good > 0){
			
					for (i=0;i<naxis1;i++){
				
						for (j=0;j<naxis2;j++){
					
							if ((radius[j][i] > 0)&&(con_dat[j][i] <= pen)&&(good_pixel(info_dat[j][i]))){
								
								/* Get the pixel location in HPL coordinates */
							
								x=rpllx+(double)i;
								y=rplly+(double)j;
								
								/* Stonyhurst/carrington needed for photospheric velocity calcs */
							
								pixel2hpl(x,y,crota2,imcrpix1,imcrpix2,cdelt1,cdelt2,&tx,&ty);	
								hpl2stony(tx,ty,dsun_obs,crln_obs,crlt_obs,&lat,&lon);
								
								/* Photospheric rotation removal */
								
								vt=v_t(lat,0,0,0)*along_phi(t,p,dsun_obs); /*0's are for default coefficients */
							
								/* Convert to spherical with no rotation to account for sat orientation and handle LOS stuff */
							
								pixel2hpl(x,y,0.0,imcrpix1,imcrpix2,cdelt1,cdelt2,&tx,&ty);
								hpl2hpc(tx,ty,&ttx,&tty);
								hpc2hcc(ttx,tty,dsun_obs,&xx,&yy,&zz);
								hcc2sphere(xx,yy,zz,&t,&p);
							
								/* LOS satellite velocity components */
							
								vr=obs_vr*along_obs(t,p,dsun_obs);
								vn=obs_vn*along_north(t,p,dsun_obs);
								vw=obs_vw*along_west(t,p,dsun_obs);
								
								/* vlos after removing all of the doppler effects. The 
								 * division by 100 is to convert cm/s to m/s. */
								
								vlos=(vlos_dat[j][i]/100)-vr-vn-vw-vt;
								
								/* Get rid of the LOS funniness */
								
								mu=sqrt(1-pow(radius[j][i]/rsun_pix,2.0));
							    vlos_fixed=vlos/limb_darken(mu);
						
								/* "Positive" sunspot polarities. */
							
								if (inc_dat[j][i] < 90){
								
									if (con_dat[j][i] <= umb){
									
										umb_p_vel[k]+=vlos;
										f_umb_p_vel[k]+=vlos_fixed;
										umb_p_counter++;
								
									}
								
									if (con_dat[j][i] > umb){
									
										f_pen_p_vel[k]+=vlos_fixed;
										pen_p_vel[k]+=vlos;
										pen_p_counter++;
									
									}
								
					
								}
					
								/* "Negative" sunspot polarities. */
					
								if (inc_dat[j][i] > 90){
							
								
									if (con_dat[j][i] <= umb){
									
										umb_n_vel[k]+=vlos;
										f_umb_n_vel[k]+=vlos_fixed;
										umb_n_counter++;
								
									}
								
									if (con_dat[j][i] > umb){
									
										f_pen_n_vel[k]+=vlos_fixed;
										pen_n_vel[k]+=vlos;
										pen_n_counter++;
									
									}
					
								}
							
							}
						}
					}
				
				}
				
				if (umb_p_counter > 0){
					
					umb_p_vel[k]/=umb_p_counter;
					f_umb_p_vel[k]/=umb_p_counter;
					
				} else {
					
					umb_p_vel[k]=0;
					f_umb_p_vel[k]=0;
					
				}
			
				if (umb_n_counter > 0){
					
					umb_n_vel[k]/=umb_n_counter;
					f_umb_n_vel[k]/=umb_n_counter;
					
				} else {
					
					umb_n_vel[k]=0;
					f_umb_n_vel[k]=0;
					
				}
				
				if (pen_p_counter > 0){
					
					pen_p_vel[k]/=pen_p_counter;
					f_pen_p_vel[k]/=pen_p_counter;
					
				} else {
					
					pen_p_vel[k]=0;
					f_pen_p_vel[k]=0;
					
				}

				if (pen_n_counter > 0){
					
					pen_n_vel[k]/=pen_n_counter;
					f_pen_n_vel[k]/=pen_n_counter;
					
				} else {
					
					pen_n_vel[k]=0;
					f_pen_n_vel[k]=0;
					
				}
				
				/* Finally, convert t_obs to a more reasonable product. */
				
				t_secs[k]=isotime(t_obs_s[k]);
			
				/* Free up dynamic arrays and close the files for the next iteration. */
				
				for (i=0;i<naxis2;i++){
				
					free(vlos_dat[i]);
					free(inc_dat[i]);
					free(con_dat[i]);
					free(radius[i]);
					free(info_dat[i]);
					
				}
				
				free(vlos_dat);
				free(inc_dat);
				free(con_dat);
				free(radius);
				free(info_dat);
				
				
				free(v_pixel_strip);
				free(i_pixel_strip);
				free(c_pixel_strip);
				free(info_pixel_strip);
			
				fits_close_file(v_ptr,&v_status);
				fits_close_file(i_ptr,&i_status);
				fits_close_file(c_ptr,&c_status);
				fits_close_file(info_ptr,&info_status);

			}
	
		}
	
		printf(" Completed run %d of %d.\n",k+1,num_records);
		
	}
	
	/* Output the results to a text file. */

	printf("Writing to file %s.\n",outfile_name);
	
	outptr=fopen(outfile_name,"w");
	bad_outptr=fopen(bad_outfile_name,"w");
	
	fprintf(outptr,"Active Region: %d\tVelocities are in m/s, fvlos is scaled to remove limb effects.\n",noaa_ar);
	fprintf(outptr,"Full UT Time            time(s)  vlos pumb  vlos ppen  vlos numb  vlos npen  fvlos pum  fvlos ppe  fvlos num  fvlos npe\n");
	fprintf(outptr,"----------------------  -------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------\n");

	fprintf(bad_outptr,"Active Region: AR%d\n",noaa_ar);
	fprintf(bad_outptr,"Full UT Time            quality\n");
	fprintf(bad_outptr,"----------------------  -------\n");
	
	for (i=0;i<num_records;i++){
		
		/* Make times relative to the first recorded one. */
		
		if (quality[i] == 0){
		
			fprintf(outptr,"%s  %7ld  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf\n",t_obs_s[i],t_secs[i]-t_secs[0],umb_p_vel[i],pen_p_vel[i],umb_n_vel[i],pen_n_vel[i],f_umb_p_vel[i],f_pen_p_vel[i],f_umb_n_vel[i],f_pen_n_vel[i]);
		
		} else {
			
			fprintf(bad_outptr,"%s  0x%08lX\n",t_obs_s[i],quality[i]);		
		
		}

	}
	
	fclose(outptr);
	fclose(bad_outptr);
	
	printf("Done!\n");
	
	/* Free memory associated with dynamic arrays. */

	free(t_secs);
	free(t_obs_s);
	
	free(umb_p_vel);
	free(pen_p_vel);
	free(umb_n_vel);
	free(pen_n_vel);
	
	free(f_umb_p_vel);
	free(f_pen_p_vel);
	free(f_umb_n_vel);
	free(f_pen_n_vel);
	
	free(quality);
	
	
	/* Close the drms records connection */
	
	drms_close_records(drms_ids,DRMS_FREE_RECORD);
    
	/* Done! */
    
	return 0;
}
