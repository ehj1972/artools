/* Calculates the mean upflow and downflow within sunspot umbrae. 
 * 
 * Changelog:		
 * 			
 * 			2014.03.23		--		First version of this code.	
 * 			2014.05.12		--		Merged in changes from other test code, handling velocity 
 * 									subtractions correctly. NOTE: We assume upflows/downflows are 
 * 									normal to the solar surface. Future code revision will put
 * 									flows along the existing magnetic field lines.						
 * 			
 */

#define _GNU_SOURCE

#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <unistd.h>
#include <dirent.h>
#include "solarlib.h"

char *module_name="sharps_upflows";
char *module_desc="determines the mean upflows/downflows within sunspot umbrae.";
char *version_id="1.0";

ModuleArgs_t module_args[]={

    {ARG_STRING,"sharps_id","","SHARPs number for the active region."},
    {ARG_STRING,"out","upflows_output.dat","name of the output file."},
	{ARG_STRING,"time_range","[]","time range to calculate flows over (optional), format is [01.02.2013_04:05:06/??]"},
	{ARG_STRING,"error","false","return error estimates on flow velocities"},
    {}

};

/* Main program */

int DoIt(){
 
    /* Variable declarations */

    /* Generic use variables */

    int i,j,k;					/* Generic loop counters */

    /* DRMS/JSOC variables */

	char *sharps_id,*time_range;				/* SHARPs ID and time range from module_args[] */
	char *outfile_name;							/* Name of the file to write out */
	char *show_error;							/* Output flags for error */
	char *drms_query;							/* The DRMS query we are going to send */
	int num_query_chars;						/* Number of characters in total query for malloc definition */
	CmdParams_t *params=&cmdparams;				/* For the module_args stuff */
	DRMS_RecordSet_t *drms_ids;					/* Holds the id structure of the records */
	DRMS_Record_t *drms_record;					/* Holds specific record information */
	DRMS_Segment_t *drms_segment;				/* Segment information */
	int drms_status;							/* Status variable */
	int num_records;							/* Number of records and variables to get the segment locations we need */
	
	/* cfitsio-related variables */

    fitsfile *d_ptr,*i_ptr,*c_ptr;		/* Pointers for FITS files */
    int d_status, i_status,c_status;	/* Status checker for opening files */
    int hdu_pos,num_hdus;				/* HDU number and location */
    int needed_hdu;						/* The HDU we want to be on */
    int any_nulls;						/* Nonzero if there are null values in data */
    char i_filename[DRMS_MAXPATHLEN+1],d_filename[DRMS_MAXPATHLEN+1];	/* Containers for the filenames */
    char c_filename[DRMS_MAXPATHLEN+1];
	

    /* Needed FITS keywords */

    int naxis1,naxis2,naxis;			/* Length of axes */
    double crpix1,crpix2;				/* SHARPs CRPIX is the center of the Sun relative to the lower left of the patch (fortran: 1,1) */
    double imcrpix1,imcrpix2;			/* SHARPs IMCRPIX is the center of the Sun in full plate coordinates (fortran) */
    double cdelt1,cdelt2;				/* HMI pixel sizes in as */
    double crota2,crln_obs,crlt_obs;	/* Rotation, carrington observer location */
    double dsun_obs,rsun_obs;			/* Distance to and radius of the sun, observed at t_obs */
    double obs_vr,obs_vw,obs_vn;		/* Satellite velocity keywords */
    double blank;						/* Values of bad data */
    TIME t_obs;							/* Observation time, ZULU */
    char t_obs_s[32];					/* String version of observation time */
    
    /* Data-related variables */
    
    double **dop_dat, **inc_dat;			/* Data arrays, size defined once naxis1/2 are known */
    double **con_dat;
    double *d_pixel_strip;					/* Strip of pixels read from the field data */
    double *i_pixel_strip;					/* Strip of pixels read from the inclination data */
    double *c_pixel_strip;					/* Strip of pixels read from the continuum data */
    long fpixel[2];							/* Holds the location of pixel to be read */
    LONGLONG nelements;						/* Number of elements in a strip to be read (naxis2) */
    FILE *outptr;							/* Pointer to the outfile */
    
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
    double lat,lon;					/* Carrington coordinate latitude and longitude */
    double p,t;						/* Spherical phi and theta coordinates */
    double x,y;						/* Pixel location on the solar disk */
    double ttx,tty,tx,ty,zz;		/* Intermediate variables for coordinate conversions */
    long *t_secs;					/* T_OBS time converted to seconds since UNIX epoch */
    double *p_flow,*n_flow;			/* Umbral flows for each polarity */
    double num_p,num_n;				/* Counters for flow means */
    double vr,vn,vw,vt;				/* Various LOS velocities due to satellite motion and solar rotation */
    
    
    /* Initialise variables that need it */

    i=j=k=0;
    drms_status=0;

	/* Grab values from module_args[] */

	sharps_id=strdup(params_get_str(params,"sharps_id"));
	outfile_name=strdup(params_get_str(params,"out"));
	time_range=strdup(params_get_str(params,"time_range"));
	show_error=strdup(params_get_str(params,"error"));

	/* Now we start the main program. The rough layout is as follows: 
	   
	   1. Query drms. Based on what we learn, either exit or grab what we need from it.
	   2. Loop through the files, calculating the centroids.
	   3. Print out a table of values at the end.				   */
	
	/* Forge the DRMS query. */
	
	num_query_chars=16+strlen(sharps_id);		/* hmi.sharp_720s[] is 16 characters */
	
	if (time_range!="[]"){
		
		num_query_chars+=strlen(time_range);
		
	}
	
	drms_query=calloc(num_query_chars,sizeof(char));
	strcat(drms_query,"hmi.sharp_720s[");
	strcat(drms_query,sharps_id);
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
	
	if (!(drms_segment_lookup(drms_record,"Dopplergram"))){
		
		printf("Dopplergram segment not present! Exiting.\n");
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
			
	t_secs=malloc(num_records*sizeof(long));
	
	p_flow=malloc(num_records*sizeof(double));
	n_flow=malloc(num_records*sizeof(double));

	for (k=0;k<num_records;k++){ 
	
		printf("Processing record %d of %d",k+1,num_records);

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
		
			/* Fill keywords */
	
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

		
			if (!(blank=drms_getkey_double(drms_record,"BLANK",&drms_status))){
					
				printf("Keyword %s not present! Exiting.\n","BLANK");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
		
			}
		
			sprint_ut(t_obs_s,t_obs);		/* Convert to string */
		
			printf(".");
		
			/* Ok, now we get the file locations for the data from DRMS. Originally I attempted
			 * to pull the data directly from DRMS but there were vague free() and malloc() issues
			 * possibly due to icc and multithreading, or library problems. */
		
			if (!(drms_segment=drms_segment_lookup(drms_record,"Dopplergram"))){
			
				printf("Problem opening the Dopplergram segment! Exiting.\n");
				drms_close_records(drms_ids,DRMS_FREE_RECORD);
				return 0;
			
			}
		
			drms_segment_filename(drms_segment,d_filename);
			
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
			
			printf(".");
		
			/* Now open the FITS files and get to work. */
	
			c_status=i_status=d_status=0;

			if (fits_open_file(&d_ptr,d_filename,READONLY,&d_status)){
				
				printf("Cannot open %s! Exiting.\n",d_filename);
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

			printf(".");

			/* Walk through the HDUs until we get to one with NAXIS=2 */

			/* Number of headers and which one we are on */

			fits_get_num_hdus(d_ptr,&num_hdus,&d_status);
			fits_get_hdu_num(d_ptr,&hdu_pos);

			/* Find the one with two naxes */

			naxis=0;
			
			for (i=1;i<=num_hdus;i++){
				
				fits_movabs_hdu(d_ptr,i,NULL,&d_status);
				fits_read_key(d_ptr,TINT,"NAXIS",&naxis,NULL,&d_status);
				
				if (naxis == 2){
	
					needed_hdu=i;
	
				}
			}

			/* Set all to the needed HDU */

			if (naxis == 0){
		
				printf("HDU problems: can't find one with the required number of axes.\n");
				fits_close_file(d_ptr,&d_status);
				fits_close_file(i_ptr,&i_status);
				fits_close_file(c_ptr,&c_status);
				exit(0);

			} else {

				fits_movabs_hdu(d_ptr,needed_hdu,NULL,&d_status);
				fits_movabs_hdu(i_ptr,needed_hdu,NULL,&i_status);
				fits_movabs_hdu(c_ptr,needed_hdu,NULL,&c_status);

			}

			printf(".");

			/* Now set some definite boundaries on arrays to help with memory usage. Arrays are called as con_dat[y][x]. */
	
			dop_dat=malloc(naxis2*sizeof(double *));
			inc_dat=malloc(naxis2*sizeof(double *));
			con_dat=malloc(naxis2*sizeof(double *));
			radius=malloc(naxis2*sizeof(double *));
			d_pixel_strip=malloc(naxis1*sizeof(double));
			i_pixel_strip=malloc(naxis1*sizeof(double));
			c_pixel_strip=malloc(naxis1*sizeof(double));
	
			for (j=0;j<naxis2;j++){
		
				dop_dat[j]=malloc(naxis1*sizeof(double));
				inc_dat[j]=malloc(naxis1*sizeof(double));
				con_dat[j]=malloc(naxis1*sizeof(double));
				radius[j]=malloc(naxis1*sizeof(double));
		
			}
    
			/* Next, put data in the arrays, clear out blank values and NaNs from the data */
	
			fpixel[0]=1;		/* cfitsio uses fortran reference instead of c when accessing data */
			nelements=naxis1;
    
			for (j=0;j<naxis2;j++){
		
				fpixel[1]=j+1;	/* Add 1 to account for fortran FITS indexing */
				fits_read_pix(d_ptr,TDOUBLE,fpixel,nelements,0,d_pixel_strip,&any_nulls,&d_status);
				fits_read_pix(i_ptr,TDOUBLE,fpixel,nelements,0,i_pixel_strip,&any_nulls,&i_status);
				fits_read_pix(c_ptr,TDOUBLE,fpixel,nelements,0,c_pixel_strip,&any_nulls,&c_status);
		
				for (i=0;i<naxis1;i++){
			
					/* Kill blank values if present, if not assign to the correct place in array */
			
					if (d_pixel_strip[i] == blank){
				
						dop_dat[j][i]=0.0;
					
					} else {
				
						dop_dat[j][i]=d_pixel_strip[i];
			
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
			
					/* Take advantage of this loop to set elements in radius to 0 */
			
					radius[j][i]=0.0;
					
				}
			}
		
			printf(".");
	
			/* We now have filled arrays of data about the SHARPs patch and all of the needed keywords. */
	
			/* Coordinate variable definitions for future use. See the declaration section */
			/* at the front of the code to understand what these refer to. */
	
			rscx=(imcrpix1-1.0);
			rscy=(imcrpix2-1.0);
			rllx=(crpix1-1.0);
			rlly=(crpix2-1.0);
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
			
			/* Determine if we actually have an active region here. */
			
			spotbounds(con_dat,naxis1,naxis2,&umb,&pen);
	
			num_good=0;
	
			for (i=0;i<naxis1;i++){
		
				for (j=0;j<naxis2;j++){
		
					if ((radius[j][i] > 0)&&(con_dat[j][i] <= pen)&&(con_dat[j][i] > 0)){
				
						num_good++;
			
					}
				}
			}
			
			/* Velocity subtractions and flow speeds */

			num_p=num_n=0;
			p_flow[k]=0;
			n_flow[k]=0;
			
			if (num_good > 0){
				
				/* We have a region. */
				
				for (i=0;i<naxis1;i++){
					
					for (j=0;j<naxis2;j++){
				
						if ((radius[j][i] > 0)&&(con_dat[j][i] <= umb)){
						
							/* There are umbrae. */
							
							/* Get the pixel location in HPL coordinates */
							
							x=rpllx+(double)i;
							y=rplly+(double)j;
							
							pixel2as(x,y,crota2,imcrpix1,imcrpix2,cdelt1,cdelt2,&tx,&ty);
							
							/* Stonyhurst/carrington needed for photospheric velocity calcs */
							
							hpl2stony(tx,ty,dsun_obs,crln_obs,crlt_obs,&lat,&lon);
							
							/* Convert to spherical with no rotation to account for sat orientation and handle LOS stuff */
							
							pixel2as(x,y,0.0,imcrpix1,imcrpix2,cdelt1,cdelt2,&tx,&ty);
							hpl2hpc(tx,ty,&ttx,&tty);
							hpc2hcc(ttx,tty,dsun_obs,&xx,&yy,&zz);
							hcc2sphere(xx,yy,zz,&t,&p);
							
							/* LOS satellite velocity components */
							
							vr=obs_vr*along_obs(t,p,dsun_obs);
							vn=obs_vn*along_north(t,p,dsun_obs);
							vw=obs_vw*along_west(t,p,dsun_obs);
							
							vt=v_t(lat)*along_phi(t,p,dsun_obs);
							
							/* Finally, subtract all of the satellite and differential rotation effects, giving a more or less 
							 * corrected data point. Divide by along_r, assuming that flows within the umbrae are normal to the
							 * solar surface.*/
							
							dop_dat[j][i]=(dop_dat[j][i]-vr-vn-vw-vt)/along_r(t,p,dsun_obs);
							
							
							/* "Positive" sunspot polarities. */
					
							if (inc_dat[j][i] > 90){
						
								p_flow[k]+=dop_dat[j][i]*is_normal(t,p,dsun_obs);
								num_p++;
					
							}
					
							/* "Negative" sunspot polarities. */
					
							if (inc_dat[j][i] < 90){
							
								n_flow[k]+=dop_dat[j][i]*is_normal(t,p,dsun_obs);
								num_n++;
					
							}
						}
					}
				}
		
			}

			if (num_p > 0){
				
				p_flow[k]/=num_p;
				
			} else {
				
				p_flow[k]=0.0;
				
			}
			
			if (num_n > 0){
				
				n_flow[k]/=num_n;
				
			} else {
				
				n_flow[k]=0.0;
				
			}

			/* Finally, convert t_obs to a more reasonable product. */
				
			t_secs[k]=isotime(t_obs_s);
			
			/* Free up dynamic arrays and close the files for the next iteration. */
				
			free(dop_dat);
			free(con_dat);
			free(inc_dat);
			free(radius);
			free(d_pixel_strip);
			free(i_pixel_strip);
			free(c_pixel_strip);
			fits_close_file(d_ptr,&d_status);
			fits_close_file(i_ptr,&i_status);
			fits_close_file(c_ptr,&c_status);
			printf(".\n");

		}
	
	}
	
	printf(".\n");
	
	/* Output the results to a text file. */

	printf("Writing to file %s.\n",outfile_name);
	
	outptr=fopen(outfile_name,"w");
	
	fprintf(outptr,"time(s)  p flow   n flow\n");
	fprintf(outptr,"-------  -------  -------\n");
	
	for (i=0;i<num_records;i++){
		
		/* Make times relative to the first recorded one. */
		
		fprintf(outptr,"%7ld  %6.2lf  %6.2lf\n",t_secs[i]-t_secs[0],p_flow[i],n_flow[i]);
	
	}
	
	fclose(outptr);
	
	printf("Done!\n");
	
	/* Free memory associated with dynamic arrays. */

	free(t_secs);
	
	free(p_flow);
	free(n_flow);
	
	/* Close the drms records connection */
	
	drms_close_records(drms_ids,DRMS_FREE_RECORD);
    
	/* Done! */
    
	return 0;
}
