/* ==================================================================
edMODISA.c

Author: Maxime Benoit-Gagne - CERC Babin and Takuvik - Canada.
Date of creation: October 3 2014.
Modified March 16 2017 by Marti Gali to use new LUT with 10 ozone levels
Usage:
./edMODISA infile lat lon doy
where
infile is the infile containing the following data set:
       Moderate Resolution Imaging Spectroradiometer (MODIS)/Aqua Aerosol 
       Cloud Water Vapor Ozone Daily L3 Global 1Deg CMG (MYD08 D3).
       This data set comes from the National Aeronautics and Space
       Administration (NASA) Level 1 and Atmosphere Archive and
       Distribution System (LAADS) (v5.1).
       The platform is Aqua.
       The sensor is MODIS.
       The temporal resolution is daily.
       The spatial resolution is 1 x 1 degree.
       Information on this data set can be found at this website:
       http://ladsweb.nascom.nasa.gov/data/ftp_site.html
       Data itself can be found at this ftp url:
       ftp://ladsweb.nascom.nasa.gov/allData/51/MYD08_D3/
lut is the lookup table that takes the wavelength, the thetas, 
    the ozone and the optical thickness to retrieve the downward
    irradiance.
    It can be a lookup table for the downward irradiance just above the
    water surface or for the downward irradiance just below the water 
    surface.
lat is the latitude in degrees N from 45 to 90.
       Return an error is lat < 45 degrees N.
       TODO: Support lat < 45 degrees N.
lon is the longitude in degrees E from -180 to 180.
doy is the day of year.

Compute the downward irradiance from a file containing atmospheric parameters
from MODIS-AQUA.
   ================================================================== */

/* ====================== INCLUDE ====================== */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "function_read_L3BIN_MODIS_Atmosphere.h"

/* ====================== SYMBOLIC CONSTANTS ====================== */
#define NBTIME 8
#define NBWL 83

/* ====================== ARRAY OF CONSTANTS ====================== */
  /*
   * Array of dimensions 8.
   * The first dimension is the index of the hours from 0 h to 21 h by step of 
   * 3 h.	       
   * The values are the hours UTC.
   * Units: h.
   */
const float ARRAY1D_ITIME_H[NBTIME] = {0., 3., 6., 9., 12., 15., 18., 21.};

/* ====================== EXTERNAL FUNCTIONS ====================== */
// In subroutine_Ed0moins_at_pixel.f
extern void ed0moins_(int *jday,
		      float *rtime,
		      float *lat,
		      float *lon,
		      float *o3,
		      float *tcl,
		      float *cf,
		      float Ed_pixel[NBWL],
		      float *ed_lut,
		      float *thetas);
extern void read_ed0moins_lut_(char *lut_fic,
			       float *ed_lut);
extern void sunpos_(int *jday,
		    float *rtime,
		    float *lat,
		    float *lon,
		    float *thetas,
		    float *azim);

/*====================== PROTOTYPES ====================== */
int get_ilat(float lat);
int get_ilon(float lon);

/* ====================== FUNCTIONS ====================== */

/*
 * IN
 * array2d_itime_iwl_ed:
 *  Array of dimensions 8 * 83.					       
 *  The first dimension is the index of the hours from 0 h to 21 h by step of 
 *  3 h. The units are h UTC.						       
 *  The second dimension is the wavelengths from 290 to 700 by step of 5 nm.  
 *  The values are the downward irradiances computed using the lookup table.
 *  Units: umol photons.m^-2.s^-1.nm^-1.                                      
 * Print on stdin.
 */
void display(float array2d_itime_iwl_ed[NBTIME][NBWL]){
  /*
   * The downward irradiances.
   * Units: umol photons.m^-2.s^-1.nm^-1.
   */
  float ed;
  // Index of the hours from 0 h to 21 h by step of 3 h. The units are h UTC.
  int itime;
  // Index of the wavelengths from 290 to 700 by step of 5 nm.
  int iwl;
  for(iwl = 0; iwl < NBWL; iwl++){
    for(itime = 0; itime < NBTIME; itime++){
      ed = array2d_itime_iwl_ed[itime][iwl];
      printf("%f ", ed);
    }
    printf("\n");
  }
}
/* ------------------------------------------------------------------ */

/*
 * IN
 * lat: The latitude in degrees N (from -90 to 90).
 * OUT
 * Return the index of the latitude on the grid of the atmospheric products
 * for MODIS-AQUA.
 */
int get_ilat(float lat){
  int ilat;
  ilat = (int)lroundf(89.5 - lat);
  return ilat;
}
/* ------------------------------------------------------------------ */

/*
 * IN
 * lon: The longitude in degrees E (from -180 to 180).
 * OUT
 * Return the index of the longitude on the grid of the atmospheric products
 * for MODIS-AQUA.
 */
int get_ilon(float lon){
  int ilon;
  ilon = (int)lroundf(179.5 + lon);
  return ilon;
}

/* ====================== MAIN ====================== */

/*
 * IN
 * argv[1]: infile. The infile containing the following data set:
 *          Moderate Resolution Imaging Spectroradiometer (MODIS)/Aqua Aerosol 
 *          Cloud Water Vapor Ozone Daily L3 Global 1Deg CMG (MYD08 D3).
 *          This data set comes from the National Aeronautics and Space
 *          Administration (NASA) Level 1 and Atmosphere Archive and
 *          Distribution System (LAADS) (v5.1).
 *          The platform is Aqua.
 *          The sensor is MODIS.
 *          The temporal resolution is daily.
 *          The spatial resolution is 1 x 1 degree.
 *          Information on this data set can be found at this website:
 *          http://ladsweb.nascom.nasa.gov/data/ftp_site.html
 *          Data itself can be found at this ftp url:
 *          ftp://ladsweb.nascom.nasa.gov/allData/51/MYD08_D3/
 * argv[2]: The lookup table that takes the wavelength, the thetas, the ozone
 *          and the optical thickness to retrieve the downward irradiance.
 *	    It can be a lookup table for the downward irradiance just above 
 *          the water surface or for the downward irradiance just below the 
 *          water surface.
 * argv[3]: lat. The latitude in degrees N from 45 to 90.
 *          Return an error is lat < 45 degrees N.
 *          TODO: Support lat < 45 degrees N.
 * argv[4]: lon. The longitude in degrees E from -180 to 180.
 * argv[5]: doy. The day of year.
 * OUT
 * Display a matirx to stdout (the Terminal).
 * Return a matrix.
 * The first dimension is the indices of the wavelengths from 290 nm to 700 nm
 * by step of 5 nm.
 * The second dimension is the indices of the times from 0h UTC (Coordinated
 * Universal Time) to 21h by step of 3 h.
 * The values are the downward irradiance computed using the lookup table.
 * The units are umol photons.m^-2.s^-1.nm^-1.
 */

int main(int argc, char *argv[]){

  /////////// Verification of the arguments. ///////////

  char* errmsg = "Usage:\n"
    "./edMODISA infile lut lat lon doy\n"
    "where\n"
    "infile is the infile containing the following data set:\n"
    "       Moderate Resolution Imaging Spectroradiometer (MODIS)/Aqua Aerosol\n"
    "       Cloud Water Vapor Ozone Daily L3 Global 1Deg CMG (MYD08 D3).\n"
    "       This data set comes from the National Aeronautics and Space\n"
    "       Administration (NASA) Level 1 and Atmosphere Archive and\n"
    "       Distribution System (LAADS) (v5.1).\n"
    "       The platform is Aqua.\n"
    "       The sensor is MODIS.\n"
    "       The temporal resolution is daily.\n"
    "       The spatial resolution is 1 x 1 degree.\n"
    "       Information on this data set can be found at this website:\n"
    "       http://ladsweb.nascom.nasa.gov/data/ftp_site.html\n"
    "       Data itself can be found at this ftp url:\n"
    "       ftp://ladsweb.nascom.nasa.gov/allData/51/MYD08_D3/\n"
    "lut is the lookup table that takes the wavelength, the thetas, \n"
    "    the ozone and the optical thickness to retrieve the downward\n"
    "    irradiance.\n"
    "    It can be a lookup table for the downward irradiance just above the\n"
    "    water surface or for the downward irradiance just below the water \n"
    "    surface.\n"
    "lat is the latitude in degrees N from 45 to 90.\n"
    "       Return an error is lat < 45 degrees N.\n"
    "lon is the longitude in degrees E from -180 to 180.\n"
    "doy is the day of year.\n";

  if(argc != 6){
    printf("%s", errmsg);
    exit(-1);
  }

  /////////// Declaration of variables. ///////////

  /*
   * Array of dimension 83.
   * The first dimension is the wavelengths from 290 to 700 by step of 5 nm.
   * The values are the downward irradiances computed using the lookup table.
   * Units: umol photons.m^-2.s^-1.nm^-1.
   */
  float array1d_iwl_ed[NBWL];
  /*
   * Array of dimensions 46 x 360.				    
   * The first dimension is the latitude from 89.5 N to 44.5 N.	    
   * The second dimension is the longitude from -179.5 E to 179.5 E.
   * The values are the cloud fraction day mean read on the 	    
   * MODIS-Atmosphere file north of 45 degrees North.               
   * No units.
   * Valid values: 0 to 1.
   * Fill value: -999.
   */
  float array2d_ilat_ilon_cf[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE];
  /*
   * Array of dimensions 46 x 360.
   * The first dimension is the latitude from 89.5 N to 44.5 N.
   * The second dimension is the longitude from -179.5 E to 179.5 E.
   * The values are the total ozone mean read on the 
   * MODIS-Atmosphere file north of 45 degrees North.
   * Units: Dobson units.
   * Fill value: -999.
   */
  float array2d_ilat_ilon_ozone[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE];
  /*
   * Array of dimensions 46 x 360.
   * The first dimension is the latitude from 89.5 N to 44.5 N.
   * The second dimension is the longitude from -179.5 E to 179.5 E.
   * The values are the cloud optical thickness combined mean read on the 
   * MODIS-Atmosphere file north of 45 degrees North.
   * No units.
   * Fill value: -999.
   */
  float array2d_ilat_ilon_taucl[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE];
  /*
   * Array of dimensions 8 * 83.					       
   * The first dimension is the index of the hours from 0 h to 21 h by step of 
   * 3 h. The units are h UTC.						       
   * The second dimension is the wavelengths from 290 to 700 by step of 5 nm.  
   * The values are the downward irradiances computed using the lookup table.
   * Units: umol photons.m^-2.s^-1.nm^-1.                                      
   * This array is used for display.                                           
   */
  float array2d_itime_iwl_ed[NBTIME][NBWL];
  /*
   * The cloud fraction.
   * No units.
   */
  float cf;
  // The day of year.
  int doy;
  /*
   * The index of the latitudes.
   */
  int ilat;
  /*
   * Index of the longitudes.
   */
  int ilon;
  /*
   * The infile containing the atmospheric data from MODIS-AQUA.
   */
  char * infile;
  // Index of the hours from 0 h to 21 h by step of 3 h. The units are h UTC.
  int itime;
  /*
   * The hour UTC.
   */
  float h;
  // The latitude in degrees N from 45 to 90.
  float lat;
  // The longitude in degrees E from -180 to 180.
  float lon;
  /*
   * The lookup table that takes the wavelength, the thetas, the ozone and     
   * the optical thickness to retrieve the downward irradiance.		       
   * It can be a lookup table for the downward irradiance just above the       
   * water surface or for the downward irradiance just below the water surface.
   */
  char * lut;
  /*
   * The ozone.
   * Units: Dobson units.
   */
  float ozone;
  /*
   * The solar azimuth angle.
   */
  float phi;
  /*
   * Array of dimensions 83 * 19 * 10 * 8.
   * The first dimension the wavelengths from 290 to 700 by step of 5 nm.
   * The second dimension is solar zenith angle. Units: Degree.
   * The third dimension the the ozone. Units: Dobson units.
   * The fourth dimension is the optical thickness. No units.
   * The values are the downward irradiances computed using the lookup table.
   * Units: umol photons.m^-2.s^-1.nm^-1.
   * The array contains the data of the lookup table.
   */
  float* ptr_array4d_iwl_ithetas_iozone_itaucl_ed;
  /*
   * The optical thickness.
   * No units.
   */
  float taucl;
  /*
   * The solar zenith angle from 0 to 90 degrees.
   */
  float thetas;

  /////////// Read the arguments. ///////////
  infile = argv[1];
  lut = argv[2];
  lat = (float)atof(argv[3]);
  lon = (float)atof(argv[4]);
  doy = atoi(argv[5]);

  /////////// Read the LUT. ///////////
  ptr_array4d_iwl_ithetas_iozone_itaucl_ed
    = (float*)malloc(sizeof(float) * NBWL * 19 * 10 * 8);
  read_ed0moins_lut_(lut,
		     ptr_array4d_iwl_ithetas_iozone_itaucl_ed);

  /////////// Read the atmospheric parameters. ///////////
  read_L3BIN_MODIS_Atmosphere(infile,
			      array2d_ilat_ilon_cf,
			      array2d_ilat_ilon_ozone,
			      array2d_ilat_ilon_taucl);
  ilat = get_ilat(lat);
  ilon = get_ilon(lon);
  cf = array2d_ilat_ilon_cf[ilat][ilon];
  ozone = array2d_ilat_ilon_ozone[ilat][ilon];
  taucl = array2d_ilat_ilon_taucl[ilat][ilon];

  /////////// Compute Ed. ///////////
  for(itime = 0; itime < NBTIME; itime++){
    h = ARRAY1D_ITIME_H[itime];
    // thetas.
    sunpos_(&doy,
	    &h,
	    &lat,
	    &lon,
	    &thetas,
	    &phi);
    // Ed.
    ed0moins_(&doy,
	      &h,
	      &lat,
	      &lon,
	      &ozone,
	      &taucl,
	      &cf,
	      array1d_iwl_ed,
	      ptr_array4d_iwl_ithetas_iozone_itaucl_ed,
	      &thetas);
    memcpy(array2d_itime_iwl_ed[itime],
	   array1d_iwl_ed,
	   sizeof(float) * NBWL);
  }

  /////////// Display. ///////////
  display(array2d_itime_iwl_ed);

  /////////// Free. ///////////
  free(ptr_array4d_iwl_ithetas_iozone_itaucl_ed);

  exit(0);
}
