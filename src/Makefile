all: sharp_mfs sharp_area sharp_field sharp_spotcount sharp_flux sharp_flux_velocity

sharp_mfs: sharp_mfs.c solarlib_analysis.c solarlib_helper.c solarlib_wcs.c
	icc -Wunused-variable -g -std=c99 -openmp -I/home/jsoc/cvs/JSOC/include -I/home/jsoc/include -L/home/jsoc/cvs/JSOC/lib/linux_x86_64 -L/home/jsoc/lib/linux-x86_64 -o $@ $+ -lm -ljsoc_main -ldrms -lcfitsio -lz -lpq -pthread -ljpeg

sharp_area: sharp_area.c solarlib_analysis.c solarlib_helper.c solarlib_wcs.c
	icc -Wunused-variable -g -std=c99 -openmp -I/home/jsoc/cvs/JSOC/include -I/home/jsoc/include -L/home/jsoc/cvs/JSOC/lib/linux_x86_64 -L/home/jsoc/lib/linux-x86_64 -o $@ $+ -lm -ljsoc_main -ldrms -lcfitsio -lz -lpq -pthread -ljpeg

sharp_field: sharp_field.c solarlib_analysis.c solarlib_helper.c solarlib_wcs.c
	icc -Wunused-variable -g -std=c99 -openmp -I/home/jsoc/cvs/JSOC/include -I/home/jsoc/include -L/home/jsoc/cvs/JSOC/lib/linux_x86_64 -L/home/jsoc/lib/linux-x86_64 -o $@ $+ -lm -ljsoc_main -ldrms -lcfitsio -lz -lpq -pthread 

sharp_spotcount: sharp_spotcount.c solarlib_analysis.c solarlib_helper.c solarlib_wcs.c
	icc -Wunused-variable -g -std=c99 -openmp -I/home/jsoc/cvs/JSOC/include -I/home/jsoc/include -L/home/jsoc/cvs/JSOC/lib/linux_x86_64 -L/home/jsoc/lib/linux-x86_64 -o $@ $+ -lm -ljsoc_main -ldrms -lcfitsio -lz -lpq -pthread 

sharp_flux: sharp_flux.c solarlib_analysis.c solarlib_helper.c solarlib_wcs.c
	icc -Wunused-variable -g -std=c99 -openmp -I/home/jsoc/cvs/JSOC/include -I/home/jsoc/include -L/home/jsoc/cvs/JSOC/lib/linux_x86_64 -L/home/jsoc/lib/linux-x86_64 -o $@ $+ -lm -ljsoc_main -ldrms -lcfitsio -lz -lpq -pthread 

sharp_flux_velocity: sharp_flux_velocity.c solarlib_analysis.c solarlib_helper.c solarlib_wcs.c
	icc -Wunused-variable -g -std=c99 -openmp -I/home/jsoc/cvs/JSOC/include -I/home/jsoc/include -I/home/ejones/include -L/home/ejones/lib -L/home/jsoc/cvs/JSOC/lib/linux_x86_64 -L/home/jsoc/lib/linux-x86_64 -o $@ $+ -lm -ljsoc_main -ldrms -lcfitsio -lz -lpq -pthread -lmpfit

clean: 
	-rm sharp_mfs sharp_area sharp_field sharp_spotcount sharp_flux sharp_flux_velocity
