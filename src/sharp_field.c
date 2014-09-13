/* This determines the total field within each polarity in an active region. 
 * 
 * All code was originally written in Python/Numpy or IDL and the ported to C, and then tested 
 * accordingly. Compile and run sharps_total_field to see the command line options.
 * 
 * 
 * Changelog:		
 * 			
 * 			2014.02.11		--		First version of this code.	
 * 
 * 			2014.03.21		--		Added in limb darkening code to help isolate umbra/penumbra better.
 * 									Added output for both polarities for umbra and penumbra. 
 * 		
 * 			2014.05.29		--		Now loads SHARPs info_map segment, using this to determine whether or not 
 * 									field data is good.
 * 
 * 			2014.06.08		--		Rejects records based upon nonzero QUALITY, outputting dates to a separate file.
 * 
 * 			2014.08.09		--		Added options to change sunspot cutoffs 		
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
	{ARG_DOUBLE,"cutoff","220","Lower limit cutoff for magnetic field in Gauss, below which is ignored."},
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
    int f_status,i_status,c_status;		/* Status checker for opening files */
    int info_status;					
    int hdu_pos,num_hdus;				/* HDU number and location */
    int needed_hdu;						/* The HDU we want to be on */
    int any_nulls;						/* Nonzero if there are null values in data */
    char i_filename[DRMS_MAXPATHLEN+1],f_filename[DRMS_MAXPATHLEN+1];		/* Containers for the filenames */
	char c_filename[DRMS_MAXPATHLEN+1],info_filename[DRMS_MAXPATHLEN+1];
	double nulval;							/* Container for what null values in data should be set to */
	long l_nulval;

    /* Needed FITS keywords */

    int naxis1,naxis2,naxis;			/* Length of axes */
    double crpix1,crpix2;				/* SHARPs CRPIX is the center of the Sun relative to the lower left of the patch (fortranL 1,1) */
    double cdelt1;						/* HMI pixel sizes in as */
    double rsun_obs;					/* Distance to and radius of the sun, observed at t_obs */
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
    FILE *outptr,*bad_outptr;				/* Pointer to the outfile */
    char outfile_name[25];					/* Name of the outfile */
    char bad_outfile_name[24];				/* Name of the outfile for bad quality*/
    
    /* Variables related to the actual task */
    
    double rllx,rlly;				/* Disk center, patch coordinates */
    double **radius;				/* Array to hold where the Sun's boundaries are, column major order */
	double umb,pen;					/* Intensity cutoffs for the umbra and penumbra */
    double num_good;				/* The number of nonzero continuum values in a sunspot region */
    double mu;						/* Angle from solar center */
    double rsun_pix,xx,yy,r;		/* Variables associated with radial boundaries */
    long *t_secs;					/* T_OBS time converted to seconds since UNIX epoch */
    double *tot_p_field,*tot_n_field;	/* Total positive and negative fields */
    double *tot_field,*tot_field_mag;	/* Total field and the sum(abs()) of the fields */
    double *umb_p_field,*umb_n_field;	/* Umbral fields */
    double *pen_p_field,*pen_n_field;	/* Penumbral fields */
	int noaa_ar;					/* Active region number */
       
    /* Initialise variables that need it */

    i=j=k=0;
    drms_status=0;
    
	/* Grab values from module_args[] */

	harpnum=strdup(params_get_str(params,"harpnum"));
	time_range=strdup(params_get_str(params,"time_range"));
	cutoff=params_get_double(params,"cutoff");
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
	tot_p_field=malloc(num_records*sizeof(double));
	tot_n_field=malloc(num_records*sizeof(double));
	umb_p_field=malloc(num_records*sizeof(double));
	pen_p_field=malloc(num_records*sizeof(double));
	umb_n_field=malloc(num_records*sizeof(double));
	pen_n_field=malloc(num_records*sizeof(double));
	tot_field=malloc(num_records*sizeof(double));
	tot_field_mag=malloc(num_records*sizeof(double));

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
			
				sprintf(outfile_name,"AR%d.field-output.dat",noaa_ar);
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
		
			if (!(rsun_obs=drms_getkey_double(drms_record,"RSUN_OBS",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","RSUN_OBS");
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
				any_nulls=0;

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
				radius=malloc(naxis2*sizeof(double *));
				info_dat=malloc(naxis2*sizeof(long *));
			
				f_pixel_strip=calloc(naxis1,sizeof(double));
				i_pixel_strip=calloc(naxis1,sizeof(double));
				c_pixel_strip=calloc(naxis1,sizeof(double));
				info_pixel_strip=calloc(naxis1,sizeof(long));
	
				for (j=0;j<naxis2;j++){
		
					field_dat[j]=calloc(naxis1,sizeof(double));
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
	
				rllx=(crpix1-1);
				rlly=(crpix2-1);
	
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
			
	
				/* Now find the total field present. */
			
				tot_p_field[k]=0;
				tot_n_field[k]=0;
				tot_field[k]=0;
				tot_field_mag[k]=0;
				umb_p_field[k]=0;
				pen_p_field[k]=0;
				umb_n_field[k]=0;
				pen_n_field[k]=0;
			
				if (num_good > 0){
			
					for (i=0;i<naxis1;i++){
				
						for (j=0;j<naxis2;j++){
					
							if ((field_dat[j][i] > cutoff)&&(radius[j][i] > 0)&&(con_dat[j][i] <= pen)&&(good_pixel(info_dat[j][i]))){
						
								/* "Positive" sunspot polarities. */
							
								if (inc_dat[j][i] < 90){
						
									tot_p_field[k]+=field_dat[j][i];
								
									if (con_dat[j][i] <= umb){
									
										umb_p_field[k]+=field_dat[j][i];
								
									}
								
									if (con_dat[j][i] > umb){
									
										pen_p_field[k]+=field_dat[j][i];
									
									}
								
									tot_field[k]+=field_dat[j][i];
					
								}
					
								/* "Negative" sunspot polarities. */
					
								if (inc_dat[j][i] > 90){
							
									tot_n_field[k]+=(-1.0*field_dat[j][i]);
								
									if (con_dat[j][i] <= umb){
									
										umb_n_field[k]+=(-1.0*field_dat[j][i]);
								
									}
								
									if (con_dat[j][i] > umb){
									
										pen_n_field[k]+=(-1.0*field_dat[j][i]);
									
									}
								
									tot_field[k]+=(-1.0*field_dat[j][i]);
					
								}
							
								/* Total field. */
							
								tot_field_mag[k]+=fabs(field_dat[j][i]);
							}
						}
					}
				
				}
			
			
				/* Free up dynamic arrays and close the files for the next iteration. */
				
				for (i=0;i<naxis2;i++){
				
					free(field_dat[i]);
					free(inc_dat[i]);
					free(con_dat[i]);
					free(radius[i]);
					free(info_dat[i]);
					
				}
				
				free(field_dat);
				free(inc_dat);
				free(con_dat);
				free(radius);
				free(info_dat);
				
				
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
	
	/* Output the results to a text file. */
	
	/* First, convert t_obs to a more reasonable product. */
	
	for (i=0;i<num_records;i++){			
	
		t_secs[i]=isotime(t_obs_s[i]);

	}

	printf("Writing to file %s.\n",outfile_name);
	
	outptr=fopen(outfile_name,"w");
	bad_outptr=fopen(bad_outfile_name,"w");
	
	fprintf(outptr,"Active Region: %d\n",noaa_ar);
	fprintf(outptr,"Full UT Time            time(s)  pos umb    pos pen    neg umb    neg pen    tot pos    tot neg    tot field  abs(tot field)\n");
	fprintf(outptr,"----------------------  -------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  --------------\n");

	fprintf(bad_outptr,"Active Region: AR%d\n",noaa_ar);
	fprintf(bad_outptr,"Full UT Time            quality  time(s)  number\n");
	fprintf(bad_outptr,"----------------------  -------  -------  ------\n");
	
	
	for (i=0;i<num_records;i++){
		
		/* Make times relative to the first recorded one. */
		
		if (quality[i] == 0){
		
			fprintf(outptr,"%s  %7ld  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf  %6.2lf\n",t_obs_s[i],t_secs[i]-t_secs[0],umb_p_field[i],pen_p_field[i],umb_n_field[i],pen_n_field[i],tot_p_field[i],tot_n_field[i],tot_field[i],tot_field_mag[i]);
		
		} else {
			
			fprintf(bad_outptr,"%s  0x%08lX  %7ld  %d\n",t_obs_s[i],quality[i],t_secs[i],i+1);	
		
		}

	}
	
	fclose(outptr);
	fclose(bad_outptr);
	
	printf("Done!\n");
	
	/* Free memory associated with dynamic arrays. */

	free(t_secs);
	free(t_obs_s);
	
	free(tot_p_field);
	free(tot_n_field);
	
	free(tot_field);
	free(tot_field_mag);
	
	free(umb_p_field);
	free(pen_p_field);
	free(umb_n_field);
	free(pen_n_field);
	
	free(quality);
	
	
	/* Close the drms records connection */
	
	drms_close_records(drms_ids,DRMS_FREE_RECORD);
    
	/* Done! */
    
	return 0;
}
