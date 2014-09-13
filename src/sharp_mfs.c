/* Footpoint separation and tilt angle code, using total field and inclination data
 * from SHARPs to determine polarity, centroid location, and footpoint separation.
 * 
 * Centroids are determined via:	
 * 					
 * 					1.	Using field and inclination data (>220 Mx default cutoff) to determine whether or not
 * 						a given data point should add to a centroid.
 * 					2.	Calculate a weighted average based on (1) to determine centroid locations. 
 * 					3.	Conversion from pixel locations to Carrington coordinates 
 * 						(see the ubiquitous Thompson, 2005 paper) for centroid locations. 
 * 									
 * The actual separation length is calculated using spherical trigonometry instead of regular
 * trigonometry. 
 * 
 * All code was originally written in Python/Numpy or IDL and the ported to C, and then tested 
 * accordingly. Compile and run sharps_mfs to see the command line options.
 * 
 * 
 * Changelog:		
 * 			
 * 			2014.02.11		--		First version of this code.
 * 
 * 			2014.05.22		--		Changed some of the output to reflect which active region we are looking at.
 * 									Outputs Carrington coordinates instead of Stonyhurst now.
 * 									Look at mean centroid latitudes, longitudes to determine hemisphere and 
 * 									the polarity of the preceeding and follower sunspots, which is output in the
 * 									header of the output file.	
 * 
 * 			2014.05.24		--		Added jpeg output options to code.
 * 			
 * 			2014.05.28		--		Array-ized T_OBS and included QUALITY keyword in output file. The QUALITY keyword
 * 									is only set when something has happened, so if it isn't present when it is pulled from
 * 									DRMS then it is set to 0.
 * 			
 * 			2014.05.29		--		Now loads SHARPs info_map segment, using this to determine whether or not 
 * 									field data is good.
 * 
 * 			2014.06.08		--		Drops everything with a QUALITY that is nonzero, outputting the times to an additional file.
 * 
 * 			2014.06.22		--		Minor bugfix involving inclinations. Note: 0-90 is positive, 90-180 negative
 * 
 * 			2014.08.08		--		Added output showing what percentage of datapoints are rejected based upon bad inversion flags
 * 
 * 			2014.08.09		--		Added options to change sunspot cutoffs
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

char *module_name="sharp_mfs";
char *module_desc="determines footpoint separation in an active region using SHARPs data. Outputs data to ARnumber.mfs-output.dat.";
char *version_id="1.0";

ModuleArgs_t module_args[]={

    {ARG_STRING,"harpnum","","SHARPs number for the active region."},
	{ARG_STRING,"time_range","[]","time range to generate footpoints over (optional), format is [01.02.2013_04:05:06/??]."},
	{ARG_DOUBLE,"cutoff","220","Lower limit cutoff for magnetic field in Gauss, below which is ignored."},
	{ARG_FLAG,"j","","When set, outputs continuum jpeg images, with p and f labels on centroid locations, with a line drawn between them. Image filename format is nnnnn.jpeg, where nnnnn is numerical ordering, starting with 00001.jpg."},
	{ARG_FLAG,"p","","When set, returns pixel locations of centroids within patch regions in the output file."},
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
	int make_jpegs,show_pixloc;					/* Output flags for jpeg creation and pixel locations */
	double cutoff;								/* Container for lower level field strength cutoff */
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

    fitsfile *f_ptr,*i_ptr,*c_ptr;		/* Pointers for FITS files */
    fitsfile *info_ptr;					
    int f_status, i_status,c_status;	/* Status checker for opening files */
    int info_status;
    int hdu_pos,num_hdus;				/* HDU number and location */
    int needed_hdu;						/* The HDU we want to be on */
    int any_nulls;						/* Nonzero if there are null values in data */
    char i_filename[DRMS_MAXPATHLEN+1],f_filename[DRMS_MAXPATHLEN+1];	/* Containers for the filenames */
    char c_filename[DRMS_MAXPATHLEN+1],info_filename[DRMS_MAXPATHLEN+1];
	double nulval;							/* Container for what null values in data should be set to */
	long l_nulval;


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
    long *quality;						/* Data quality keyword */
    
    /* Data-related variables */
    
    double **field_dat, **inc_dat;			/* Data arrays, size defined once naxis1/2 are known */
    double **con_dat;
    long **info_dat;
    double *f_pixel_strip;					/* Strip of pixels read from the field data */
    double *i_pixel_strip;					/* Strip of pixels read from the inclination data */
    double *c_pixel_strip;					/* Strip of pixels read from the continuum data */
    long *info_pixel_strip;					/* Strip of pixels read from the continuum data */
    long fpixel[2];							/* Holds the location of pixel to be read */
    LONGLONG nelements;						/* Number of elements in a strip to be read (naxis2) */
    FILE *outptr,*bad_outptr;							/* Pointers to the data and quality outfiles */
    char outfile_name[24];					/* Name of the outfile for data */
    char bad_outfile_name[24];				/* Name of the outfile for bad quality*/
    
    /* Variables related to the actual task */
    
    double rscx,rscy;				/* Disk center, plate coordinates */
    double rllx,rlly;				/* Disk center, patch coordinates */
    double rpllx,rplly;				/* Lower left of patch, plate coordinates */
    double umb,pen;					/* Intensity cutoffs for the umbra and penumbra */
    double num_good;				/* The number of nonzero continuum values in a sunspot region */
    double max_con_value;			/* Maximum continuum intensity value in patch */
    double mu;						/* Angle from solar center */
    double **radius;				/* Array to hold where the Sun's boundaries are, column major order */
    double rsun_pix,xx,yy,r;		/* Variables associated with radial boundaries */
    double *pcx,*pcy;				/* Pixel locations of the positive centroid */
    double *ncx,*ncy;				/* Pixel locations of the negative centroid */
    double p_tot,n_tot;				/* Total positive and negative pixels in region */
    double p_x,p_y,n_x,n_y;			/* Centroid locations in plate coordinates */
    double *p_lat,*p_lon;			/* Carrington coordinates for positive centroid */
    double *n_lat,*n_lon;			/* Carrington coordinates for negative centroid */
    double ptx,pty,ntx,nty;			/* Intermediate variables for coordinate conversions */
    long *t_secs;					/* T_OBS time converted to seconds since UNIX epoch */
    double *length,*tilt;			/* Footpoint separation and tilt angle */
    double weight;					/* Placeholder variable for field strenth */
	double preceding;				/* Counter to help determine which polarity is preceding/follower */
	double hemisphere;				/* Counter for hemisphere */
	int noaa_ar;					/* Active region number */
	char leading_spots[10];			/* Character array for keeing track of preceding/follower spots */
    char which_hemisphere[6];		/* Keeps track of which hemisphere this is */
    double pos_bad,neg_bad;			/* Counters for tracking the number of bad positive/negative flags in field data */
    double pos_good,neg_good;		/* Track good ones */
    
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
    preceding=0;
    hemisphere=0;

	/* Grab values from module_args[] */

	harpnum=strdup(params_get_str(params,"harpnum"));
	time_range=strdup(params_get_str(params,"time_range"));
	make_jpegs=params_isflagset(params,"j");
	show_pixloc=params_isflagset(params,"p");
	cutoff=params_get_double(params,"cutoff");
	umb_sig=params_get_double(params,"umb_sig");
	pen_sig=params_get_double(params,"pen_sig");
	
	/* Now we start the main program. The rough layout is as follows: 
	   
	   1. Query drms. Based on what we learn, either exit or grab what we need from it.
	   2. Loop through the files, calculating the centroids.
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
	

	if (!(drms_segment_lookup(drms_record,"field"))){
		
		printf("Field segment not present! Exiting.\n");
		drms_close_records(drms_ids,DRMS_FREE_RECORD);
		return 0;
	
	}

	if (!(drms_segment_lookup(drms_record,"inclination"))){
		
		printf("Inclination segment not present! Exiting.\n");
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
	
	/* Pixel locations */
	
	pcx=malloc(num_records*sizeof(double));
	pcy=malloc(num_records*sizeof(double));
	ncx=malloc(num_records*sizeof(double));
	ncy=malloc(num_records*sizeof(double));
	
	/* Carrington/Stonyhurst locations */
	
	p_lat=malloc(num_records*sizeof(double));
	p_lon=malloc(num_records*sizeof(double));
	n_lat=malloc(num_records*sizeof(double));
	n_lon=malloc(num_records*sizeof(double));
	
	/* Time-related quantities */
	
	t_secs=malloc(num_records*sizeof(long));
	
	t_obs_s=malloc(num_records*sizeof(char *));
	
	for (i=0;i<num_records;i++){
		
		t_obs_s[i]=malloc(T_OBS_LENGTH*sizeof(char));
		
	}
		
	/* Calculation results */
	
	length=malloc(num_records*sizeof(double));
	tilt=malloc(num_records*sizeof(double));
	
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
				
				/* Create data and bad record output filenames */
			
				sprintf(outfile_name,"AR%d.mfs-output.dat",noaa_ar);
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
			
				if (make_jpegs){
					
					printf("We are creating images.\n");
				
				} else {
			
					printf("We are not creating images.\n");
				
				}

				if (show_pixloc){
					
					printf("Dataset output will contain pixel location of centroids.\n");
				
				} else {
			
					printf("Default dataset output.\n");
				
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
			
			if (quality[k] == 0){
		
				/* Ok, now we get the file locations for the data from DRMS. Originally I attempted
				* to pull the data directly from DRMS but there were vague free() and malloc() issues
				* possibly due to icc and multithreading, or library problems. */
		
				if (!(drms_segment=drms_segment_lookup(drms_record,"field"))){
			
					printf("Problem opening the field segment! Exiting.\n");
					drms_close_records(drms_ids,DRMS_FREE_RECORD);
					return 0;
			
				}
		
				drms_segment_filename(drms_segment,f_filename);
			
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
	
				f_status=0;
				i_status=0;
				c_status=0;
				info_status=0;

				if (fits_open_file(&f_ptr,f_filename,READONLY,&f_status)){
				
					printf("Cannot open %s! Exiting.\n",f_filename);
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

				fits_get_num_hdus(f_ptr,&num_hdus,&f_status);
				fits_get_hdu_num(f_ptr,&hdu_pos);

				/* Find the one with two naxes */

				naxis=0;
			
				for (i=1;i<=num_hdus;i++){
				
					fits_movabs_hdu(f_ptr,i,NULL,&f_status);
					fits_read_key(f_ptr,TINT,"NAXIS",&naxis,NULL,&f_status);
				
					if (naxis == 2){
	
						needed_hdu=i;
	
					}
				}

				/* Set all to the needed HDU */

				if (naxis == 0){
		
					printf("HDU problems: can't find one with the required number of axes.\n");
					fits_close_file(f_ptr,&f_status);
					fits_close_file(i_ptr,&i_status);
					fits_close_file(c_ptr,&c_status);
					fits_close_file(info_ptr,&info_status);
					exit(0);

				} else {

					fits_movabs_hdu(f_ptr,needed_hdu,NULL,&f_status);
					fits_movabs_hdu(i_ptr,needed_hdu,NULL,&i_status);
					fits_movabs_hdu(c_ptr,needed_hdu,NULL,&c_status);
					fits_movabs_hdu(info_ptr,needed_hdu,NULL,&info_status);

				}

				/* Now set some definite boundaries on arrays to help with memory usage. Arrays are called as con_dat[y][x]. */
	
				field_dat=malloc(naxis2*sizeof(double *));
				inc_dat=malloc(naxis2*sizeof(double *));
				con_dat=malloc(naxis2*sizeof(double *));
				info_dat=malloc(naxis2*sizeof(long *));
			
			
				if (make_jpegs){
			
					jpeg_outdata=malloc(naxis2*sizeof(int *));
			
				}
			
				radius=malloc(naxis2*sizeof(double *));
				
				f_pixel_strip=calloc(naxis1,sizeof(double));
				i_pixel_strip=calloc(naxis1,sizeof(double));
				c_pixel_strip=calloc(naxis1,sizeof(double));
				info_pixel_strip=calloc(naxis1,sizeof(long));
				
				row_pointer[0]=(unsigned char*)malloc(naxis1*image_bytes);
	
				for (j=0;j<naxis2;j++){
		
					field_dat[j]=calloc(naxis1,sizeof(double));
					inc_dat[j]=calloc(naxis1,sizeof(double));
					con_dat[j]=calloc(naxis1,sizeof(double));
					info_dat[j]=calloc(naxis1,sizeof(long));
					radius[j]=calloc(naxis1,sizeof(double));
				
					if (make_jpegs){
					
						jpeg_outdata[j]=calloc(naxis1,sizeof(int));
				
					}
		
				}
    
				/* Next, put data in the arrays, clear out blank values and NaNs from the data */
	
				fpixel[0]=1;		/* cfitsio uses fortran reference instead of c when accessing data */
				nelements=naxis1;
				nulval=0;
				l_nulval=0;
				any_nulls=0;
    
				for (j=0;j<naxis2;j++){
		
					fpixel[1]=j+1;	/* Add 1 to account for fortran FITS indexing */
				
					fits_read_pix(f_ptr,TDOUBLE,fpixel,nelements,&nulval,f_pixel_strip,&any_nulls,&f_status);
					fits_read_pix(i_ptr,TDOUBLE,fpixel,nelements,&nulval,i_pixel_strip,&any_nulls,&i_status);
					fits_read_pix(c_ptr,TDOUBLE,fpixel,nelements,&nulval,c_pixel_strip,&any_nulls,&c_status);
					fits_read_pix(info_ptr,TLONG,fpixel,nelements,&l_nulval,info_pixel_strip,&any_nulls,&info_status);
				
					for (i=0;i<naxis1;i++){
			
						/* Kill blank values if present, if not assign to the correct place in array */
				
						if (f_pixel_strip[i] == blank){
				
							field_dat[j][i]=0.0;
						
						} else {
				
							field_dat[j][i]=f_pixel_strip[i];
			
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
			
				if (make_jpegs){
				
					/* Write integer values with 255 being maximum for greyscale image output */
				
					for (j=0;j<naxis2;j++){
						for (i=0;i<naxis1;i++){
							jpeg_outdata[j][i]=(int)(255*con_dat[j][i]/max_con_value);
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
			
				/* Now find the pixel location of the centroids. Locations are set to -1000 if polarities aren't present. */
	
				pcx[k]=pcy[k]=0.0;
				ncx[k]=ncy[k]=0.0;
				p_tot=n_tot=0.0;
				pos_bad=neg_bad=0.0;
				pos_good=neg_good=0.0;
			
				if (num_good > 0){
				
					/* We have a region. */
				
					for (i=0;i<naxis1;i++){
					
						for (j=0;j<naxis2;j++){
				
							if ((radius[j][i] > 0)&&(field_dat[j][i] > cutoff)&&(con_dat[j][i] <= pen)){
							
								if (good_pixel(info_dat[j][i])){
						
									/* "Positive" sunspot polarities. The default 220 limit is based upon SHARPs field data 
									* below 220 Mx/cm^2 being considered noise. Checking against bad pixel data as well. */
							
									weight=field_dat[j][i];
					
									if (inc_dat[j][i] < 90){
						
										pcx[k]+=weight*(double)i;
										pcy[k]+=weight*(double)j;
										p_tot+=weight;
										pos_good++;
					
									}
					
									/* "Negative" sunspot polarities. */
					
									if (inc_dat[j][i] > 90){
							
										ncx[k]+=weight*(double)i;
										ncy[k]+=weight*(double)j;
										n_tot+=weight;
										neg_good++;
					
									}
								
								} else {
									
									if (inc_dat[j][i] < 90){
										
										pos_bad++;
										
									}
									
									if (inc_dat[j][i] > 90){
										
										neg_bad++;
										
									}
									
								}
								
							}
						}
					}
					
					/* Fractional error reporting */
					
					if ((pos_good > 0)&&(neg_good > 0)){
					
						printf(" Positive bad to good pixel percentage: %6.2f\n Negative bad to good pixel percentage: %6.2f\n",100*(pos_bad/pos_good),100*(neg_bad/neg_good));
				
					} else {
						
						printf(" No good pixels present in one centroid\n");
						
					}
				
					/* Error handling and final pixel centroid calculations, return -1000 for missing polarites. */
		
					if (p_tot == 0.0){
			
						pcx[k]=pcy[k]=-10000;
		
					} else {
			
						pcx[k]/=p_tot;
						pcy[k]/=p_tot;
		
					}
	
					if (n_tot == 0.0){
			
						ncx[k]=ncy[k]=-10000;
			
					} else {
			
						ncx[k]/=n_tot;
						ncy[k]/=n_tot;
		
					}
		
				} else {
		
					/* We have no region. */
		
					pcx[k]=pcy[k]=-10000;
					ncx[k]=ncy[k]=-10000;
	
				}
				
				
				/* Centroid locations in plate coordinates. */
	
				p_x=rpllx+pcx[k];
				p_y=rplly+pcy[k];
				n_x=rpllx+ncx[k];
				n_y=rplly+ncy[k];
	
				/* Convert plate coordinates to Carrington coordinates. */
			
				if ((p_x > 0)&&(p_y > 0)){
		
					pixel2hpl(p_x,p_y,crota2,imcrpix1,imcrpix2,cdelt1,cdelt2,&ptx,&pty);
					hpl2stony(ptx,pty,dsun_obs,crln_obs,crlt_obs,&p_lat[k],&p_lon[k]);
					
					if (p_lon[k] != -1000){		
					
						p_lon[k]-=crln_obs;		/* Stony to Carrington off by crln_obs */
						
					}
					
	
				} else {
		
					p_lat[k]=p_lon[k]=-1000;
	
				}
	
	
				if ((n_x > 0)&&(n_y > 0)){
		
					pixel2hpl(n_x,n_y,crota2,imcrpix1,imcrpix2,cdelt1,cdelt2,&ntx,&nty);
					hpl2stony(ntx,nty,dsun_obs,crln_obs,crlt_obs,&n_lat[k],&n_lon[k]);

					if (n_lon[k] != -1000){

						n_lon[k]-=crln_obs;		/* Stony to Carrington off by crln_obs */
						
					}
		
				} else {
		
					n_lat[k]=n_lon[k]=-1000;
	
				}
				
				if ((p_lon[k] > -360)&&(n_lon[k] > -360)){
				
					if (p_lon[k] > n_lon[k]){
					
						preceding+=1;
					
					} else {
					
						preceding-=1;
					}
						
					if ((0.5*(p_lat[k]+n_lat[k])) > 0){
					
						hemisphere+=1;
					
					} else {
					
						hemisphere-=1;
					}
				
				}
			
				/* Calculate footpoint separation and tilt angles. */
			
				arclength(p_lon[k],n_lon[k],p_lat[k],n_lat[k],&length[k]);
			
				tilt_angle(p_lon[k],n_lon[k],p_lat[k],n_lat[k],&tilt[k]);
			
				/* Write JPEG */
			
				if (make_jpegs){
					
					printf(" Creating jpeg image.\n");
				
					/* Need to draw circles where the centroids are, and then a line between them */
				
					if ((p_x > 0)&&(p_y> 0)){
					
						draw_centroid(pcx[k],pcy[k],3.0,naxis1,naxis2,jpeg_outdata,0);
					
					}

					if ((n_x > 0)&&(n_y > 0)){
					
						draw_centroid(ncx[k],ncy[k],3.0,naxis1,naxis2,jpeg_outdata,0);
					
					}
				
					if ((n_x > 0)&&(n_y > 0)&&(p_x > 0)&&(p_y > 0)){
					
						if (ncx[k] > pcx[k]){
					
							draw_line(pcx[k],ncx[k],pcy[k],ncy[k],naxis1,naxis2,jpeg_outdata,0,2);
						
						} else {
						
							draw_line(ncx[k],pcx[k],ncy[k],pcy[k],naxis1,naxis2,jpeg_outdata,0,2);
						
						}
					}
				
					sprintf(jpg_outfile,"%05d.jpg",k+1);
			
					jpeg_ptr=fopen(jpg_outfile,"wb");
	
					cinfo.err=jpeg_std_error(&jerr);
					jpeg_create_compress(&cinfo);
					jpeg_stdio_dest(&cinfo,jpeg_ptr);
	
					cinfo.image_width=naxis1;
					cinfo.image_height=naxis2;
					cinfo.input_components=image_bytes;
					cinfo.in_color_space=colorscheme;
	
					jpeg_set_defaults(&cinfo);
					jpeg_set_quality(&cinfo,80,TRUE);
					jpeg_start_compress(&cinfo,TRUE);
	
					for (j=0;j<naxis2;j++){
						for (i=0;i<naxis1;i++){
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
				
					free(field_dat[i]);
					free(con_dat[i]);
					free(inc_dat[i]);
					free(info_dat[i]);
					free(radius[i]);
				
					if (make_jpegs){
					
						free(jpeg_outdata[i]);
				
					}
				
				}
				
				free(field_dat);
				free(con_dat);
				free(inc_dat);
				free(info_dat);
			
				if (make_jpegs){
					
					free(jpeg_outdata);
					free(row_pointer[0]);
				
				}
			
				free(f_pixel_strip);
				free(i_pixel_strip);
				free(c_pixel_strip);
				free(info_pixel_strip);
			
				fits_close_file(f_ptr,&f_status);
				fits_close_file(i_ptr,&i_status);
				fits_close_file(c_ptr,&c_status);
				fits_close_file(info_ptr,&info_status);

			}
	
		}
		printf(" Completed run %d of %d.\n",k+1,num_records);
	
	}
	
	if (preceding > 0){
		
		strcpy(leading_spots,"positive");
		
	} else {
		
		strcpy(leading_spots,"negative");
	
	}
	
	if (hemisphere > 0){
		
		strcpy(which_hemisphere,"north");
		
	} else {
		
		strcpy(which_hemisphere,"south");
	
	}	

	
	/* Output the results to a text file. */
	
	/* First, convert t_obs to a more reasonable product. */
	
	for (i=0;i<num_records;i++){			
	
		t_secs[i]=isotime(t_obs_s[i]);

	}
	
	printf("Writing to file %s.\n",outfile_name);
	
	outptr=fopen(outfile_name,"w");
	bad_outptr=fopen(bad_outfile_name,"w");
	
	if (show_pixloc){
	
		fprintf(outptr,"Active Region: AR%d\tHemisphere: %s\tPreceding polarity: %s\n",noaa_ar,which_hemisphere,leading_spots);
		fprintf(outptr,"Full UT Time            img num  time(s)  pos x    pos y    neg x    neg y    pos lat  pos lon  neg lat  neg lon  sep(km)  tilt(deg)\n");
		fprintf(outptr,"----------------------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  -------  ---------\n");
		
		fprintf(bad_outptr,"Active Region: AR%d\n",noaa_ar);
		fprintf(bad_outptr,"Full UT Time            quality  time(s)  number\n");
		fprintf(bad_outptr,"----------------------  -------  -------  ------\n");
	
	
	
		for (i=0;i<num_records;i++){
		
			/* Make times relative to the first recorded one. */
			
			if(quality[i] == 0){
		
				fprintf(outptr,"%s  %7d  %7ld  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf\n",t_obs_s[i],i+1,t_secs[i]-t_secs[0],pcx[i],pcy[i],ncx[i],ncy[i],p_lat[i],p_lon[i],n_lat[i],n_lon[i],length[i]/(1000.0),tilt[i]);
			
			} else {
				
				fprintf(bad_outptr,"%s  0x%08lX  %7ld  %d\n",t_obs_s[i],quality[i],t_secs[i],i+1);
				
			}		
				
		}
		
	} else {
		
		fprintf(outptr,"Active Region: AR%d\tHemisphere: %s\tPreceding polarity: %s\n",noaa_ar,which_hemisphere,leading_spots);
		fprintf(outptr,"Full UT Time            img num  time(s)  pos lat  pos lon  neg lat  neg lon  sep(km)  tilt(deg)\n");
		fprintf(outptr,"----------------------  -------  -------  -------  -------  -------  -------  -------  ---------\n");

	fprintf(bad_outptr,"Active Region: AR%d\n",noaa_ar);
	fprintf(bad_outptr,"Full UT Time            quality  time(s)  number\n");
	fprintf(bad_outptr,"----------------------  -------  -------  ------\n");
	
		for (i=0;i<num_records;i++){
		
			/* Make times relative to the first recorded one. */
			
			if(quality[i] == 0){
		
				fprintf(outptr,"%s  %7d  %7ld  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf\n",t_obs_s[i],i+1,t_secs[i]-t_secs[0],p_lat[i],p_lon[i],n_lat[i],n_lon[i],length[i]/(1000.0),tilt[i]);
			
			} else {
				
				fprintf(bad_outptr,"%s  0x%08lX  %7ld  %d\n",t_obs_s[i],quality[i],t_secs[i],i+1);
				
			}	
			
		}
	}
	
	
	fclose(outptr);
	fclose(bad_outptr);
	
	printf("Done!\n");
	
	/* Free memory associated with dynamic arrays. */
    

	free(pcx);
	free(pcy);
	free(ncx);
	free(ncy);
	
	free(quality);
	
	free(p_lat);
	free(p_lon);
	free(n_lat);
	free(n_lon);

	free(t_secs);
	free(t_obs_s);
	
	free(length);
	free(tilt);
	
	/* Close the drms records connection */
	
	drms_close_records(drms_ids,DRMS_FREE_RECORD);
    
	/* Done! */
    
	return 0;
}
