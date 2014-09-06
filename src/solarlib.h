/* Header file containing prototypes from solarlib.c */

/* Analysis-related prototypes */

double limb_darken(double);		/* Scaling factor to remove limb darkening effects */
void spotbounds(double**,int,int,double*,double*,double,double); /* Determines the cutoff intensity values for the umbra and penumbra */
void arclength(double,double,double,double,double*); /* Calculates the arclength between two centroids */
void arclength_e(double,double,double,double,double,double,double,double,double*,double*); /* Arclength, with error */
void tilt_angle(double,double,double,double,double*); /* Calculates the tilt angle between two centroids relative to the solar equator */
void tilt_angle_e(double,double,double,double,double,double,double,double,double*,double*); /* tilt_angle, with error */
double v_t(double);		/* Tangential velocity (along the latitude direction) from photospheric data */
double centroid_error(double*,double*,int,double,double); /* Performs centroid error calculations based upon field error */
double br_error(double,double,double,double); /* Performs brightness ratio error calculations based on continuum intensity information */
int connected_component_detect(int **,int,int); /* Finds the number of unique sunspots within an image mask */
	
/* Helper prototypes */

double arg(double,double);		/* The arg(x,y) function as defined in Thompson's 2005 paper */
void rotate(double,double,double,double*,double*);	/* Rotates a vector CCW by default, as defined by CROTA2 in HMI's keywords */
long isotime(char*);	/* Converts HMI zulu time ISO format to seconds for easier use */
double omega(double,double,double,double); /* Helper function for v_t */
int check_labels(int, int, int**, int, int); /* Needed by connected_component_detect */
void draw_centroid(double,double,double,int,int,int **,int); /* Draws filled circles for jpeg output */
void draw_line(double,double,double,double,int,int,int **,int,int); /* Draws line between centroid to show separation */
void lint_to_binary(long,int *);	/* Converts a long integer to binary */
int good_pixel(long);	/* Checks an error code in long int against SHARPs criteria to see if it falls within the range of good pixels */
double bilinear_interpolate(double,double,int,int,double**); /* 2-D bicubic interpolation */
double pixel_area(double,double,double);	/* Returns the area of a pixel bounded by latitude/longitude in degrees */

/* WCS-related prototypes, see Thompson's 2005 paper on solar coordinates */

void pixel2hpl(double,double,double,double,double,double,double,double*,double*);	/* Converts pixel locations to HPL-TAN/HPLN-TAN */
void hpl2pixel(double,double,double,double,double,double,double,double*,double*);	/* Converts HPL-TAN/HPLN-TAN to pixel locations */
void hpl2hpc(double,double,double*,double*);	/* Converts HPL-TAN/HPLN-TAN into helioprojective cartesian */
void hpc2hpl(double,double,double*,double*);	/* Converts helioprojective cartesian into HPL-TAN/HPLN-TAN */
void hpc2hcc(double,double,double,double*,double*,double*);	/* Converts helioprojective cartesian to heliocentric cartesian */
void hcc2hpc(double,double,double,double,double*,double*);	/* Converts heliocentric cartesian to helioprojective cartesian */
void hcc2stony(double,double,double,double,double,double*,double*);		/* Converts heliocentric cartesian to Stonyhurst */
void stony2hcc(double,double,double,double,double*,double*,double*);	/* Converts Stonyhurst to heliocentric cartesian */
void hpl2stony(double,double,double,double,double,double*,double*);		/* Converts HPL-TAN/HPLN-TAN into Stonyhurst */
void stony2hpl(double,double,double,double,double,double*,double*);		/* Converts Stonyhurst into HPL-TAN/HPLN-TAN */
void hcc2sphere(double,double,double,double*,double*);	/* Converts heliocentric cartesian to standard (physics) spherical */
double along_obs(double,double,double);	
double along_north(double,double,double);
double along_west(double,double,double);
double along_phi(double,double,double);
double along_theta(double,double,double);
double along_r(double,double,double);
double is_normal(double,double,double);

/* Cartography (located in WCS) */

void CEA2latlon(double,double,double,double*,double*);
