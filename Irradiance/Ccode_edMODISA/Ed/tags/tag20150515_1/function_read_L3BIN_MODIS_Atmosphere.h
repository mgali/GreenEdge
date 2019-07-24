/*
 * function_read_L3BIN_MODIS_Atmosphere.h
 * Maxime Benoit-Gagne
 * November 28th, 2012
 *
 * function_read_L3BIN_MODIS_Atmosphere.h is the header for the function in
 * function_read_L3BIN_MODIS_Atmosphere.c.
 */

#ifndef _L3BIN_MODIS_ATMOSPHERE_H
#define _L3BIN_MODIS_ATMOSPHERE_H
#define NBLAT_MODIS_ATMOSPHERE 46
#define NBLON_MODIS_ATMOSPHERE 360

/*
 * filename: The filename of a L3 daily file of MODIS Atmosphere.
 * CFDM: Cloud_Fraction_Day_Mean. No units. Valid values: 0 to 1. 
 *       Fill value: -999.
 * TO3M: Total_Ozone_Mean. Dobson units. Fill value: -999.
 * COTCM: Cloud_Optical_Thickness_Combined_Mean. No units. Fill value: -999.
 * Read a file in the three arrays.
 */
void read_L3BIN_MODIS_Atmosphere(char filename[],
				 float CFDM[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE],
				 float TO3M[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE],
				 float COTCM[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE]);

#endif
