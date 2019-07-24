/*
 * See comments in the header function_read_L3BIN_MODIS_Atmosphere.h.
 */

///////////////////////////////// INCLUDES /////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <mfhdf.h>

#include "function_read_L3BIN_MODIS_Atmosphere.h"

///////////////////////////////// CONSTANTS /////////////////////////////////

#define ATTR_INDEX_SCALE_FACTOR 4
#define FILL_VALUE_MODIS_ATMOSPHERE -9999
#define FILL_VALUE_TAKUVIK -999.
#define MAX_VAR_DIMS 2
#define SDS_INDEX_CLOUD_FRACTION_DAY_MEAN 347
#define SDS_INDEX_CLOUD_OPTICAL_THICKNESS_COMBINED_MEAN 417
#define SDS_INDEX_TOTAL_OZONE_MEAN 631

///////////////////////////////// FUNCTION /////////////////////////////////

void read_L3BIN_MODIS_Atmosphere(char filename[],
				 float CFDM[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE],
				 float TO3M[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE],
				 float COTCM[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE]){
  int i,j;
  int ionedata;
  float fonedata;
  intn status;                        /* function status flag */
  int32 Hid;                          /* HDF interface ID */
  int32 SDid;                         /* Scientific Data interface ID */
  int32 sds_id;                         /* Annotation interface ID */
  int32 rank, num_type,  dim_sizes[MAX_VAR_DIMS], n_attrs;
  char name[1024];
  double scale_buf;
  int16 data[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE];
  int32 start[2] = {0, 0};
  int32 edge[2]={NBLAT_MODIS_ATMOSPHERE, NBLON_MODIS_ATMOSPHERE};
  if ((Hid = Hopen(filename, DFACC_READ, 0)) == FAIL){
    fprintf(stderr, "Error: unable to open file %s\n", filename);
    exit(-1);
  }
  SDid = SDstart(filename, DFACC_READ);       /* initiate SD API interface */

  // Cloud_Fraction_Day_Mean
  sds_id = SDselect(SDid, SDS_INDEX_CLOUD_FRACTION_DAY_MEAN);
  status = SDgetinfo(sds_id, name, &rank, dim_sizes, &num_type, &n_attrs);
  status=SDreadattr(sds_id, ATTR_INDEX_SCALE_FACTOR, &scale_buf);
  //printf("%s %f \n",name, scale_buf);
  status = SDreaddata(sds_id,start,NULL,edge,(VOIDP) data);	
  for (i = 0; i < NBLAT_MODIS_ATMOSPHERE; i++){
    for (j = 0; j < NBLON_MODIS_ATMOSPHERE; j++){
      ionedata = (int)data[i][j];
      if(ionedata == FILL_VALUE_MODIS_ATMOSPHERE){
	fonedata = FILL_VALUE_TAKUVIK;
      }else{
	fonedata = (float)(ionedata * scale_buf);
      }
      CFDM[i][j] = fonedata;
    }
  }

  // Total_Ozone_Mean
  sds_id = SDselect(SDid, SDS_INDEX_TOTAL_OZONE_MEAN);
  status = SDgetinfo(sds_id, name, &rank, dim_sizes, &num_type, &n_attrs);
  status=SDreadattr(sds_id, ATTR_INDEX_SCALE_FACTOR, &scale_buf);
  //printf("%s %f \n",name, scale_buf);
  status = SDreaddata(sds_id,start,NULL,edge,(VOIDP) data);	
  for (i = 0; i < NBLAT_MODIS_ATMOSPHERE; i++){
    for (j = 0; j < NBLON_MODIS_ATMOSPHERE; j++){
      ionedata = (int)data[i][j];
      if(ionedata == FILL_VALUE_MODIS_ATMOSPHERE){
	fonedata = FILL_VALUE_TAKUVIK;
      }else{
	fonedata = (float)(ionedata * scale_buf);
      }
      TO3M[i][j] = fonedata;
    }
  }

  // Cloud_Optical_Thickness_Combined_Mean
  sds_id = SDselect(SDid, SDS_INDEX_CLOUD_OPTICAL_THICKNESS_COMBINED_MEAN);
  status = SDgetinfo(sds_id, name, &rank, dim_sizes, &num_type, &n_attrs);
  status=SDreadattr(sds_id, ATTR_INDEX_SCALE_FACTOR, &scale_buf);
  //printf("%s %f \n",name, scale_buf);
  status = SDreaddata(sds_id,start,NULL,edge,(VOIDP) data);	
  for (i = 0; i < NBLAT_MODIS_ATMOSPHERE; i++){
    for (j = 0; j < NBLON_MODIS_ATMOSPHERE; j++){
      ionedata = (int)data[i][j];
      if(ionedata == FILL_VALUE_MODIS_ATMOSPHERE){
	fonedata = FILL_VALUE_TAKUVIK;
      }else{
	fonedata = (float)(ionedata * scale_buf);
      }
      COTCM[i][j] = fonedata;
    }
  }

  status=SDendaccess(sds_id);

}
