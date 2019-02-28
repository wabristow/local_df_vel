/* make_pgrid.c
   ============ 
   Author: R.J.Barnes
*/

/*
 LICENSE AND DISCLAIMER
 
 Copyright (c) 2012 The Johns Hopkins University/Applied Physics Laboratory
 
 This file is part of the Radar Software Toolkit (RST).
 
 RST is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 any later version.
 
 RST is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with RST.  If not, see <http://www.gnu.org/licenses/>.
 
 
 
*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aacgm.h"

#include "local_df_vel.h"



int make_local_pgrid(struct cell *gptr, int ngrid, struct CnvGrid *ptr) {
  
  int i,j;

  int num;
  int poly; 

  double lat, lon;
  double mlat, mlon,r;

  ptr->type=0;
  ptr->num=0;
 
  num=0;  
  poly=0;

  int nlon=(int) (360.0/2); 

  if (ptr->vertex !=NULL) free(ptr->vertex);
  if (ptr->lat !=NULL) free(ptr->lat);
  if (ptr->lon !=NULL) free(ptr->lon);
  if (ptr->mag !=NULL) free(ptr->mag);
  if (ptr->azm !=NULL) free(ptr->azm);
  if (ptr->ex !=NULL) free(ptr->ex);
  if (ptr->ey !=NULL) free(ptr->ey);

  ptr->vertex=malloc(sizeof(int)*ngrid*4);
  ptr->lat=malloc(sizeof(double)*ngrid);
  ptr->lon=malloc(sizeof(double)*ngrid);
  ptr->mag=malloc(sizeof(double)*ngrid);
  ptr->azm=malloc(sizeof(double)*ngrid);
  ptr->ex=malloc(sizeof(double)*ngrid);
  ptr->ey=malloc(sizeof(double)*ngrid);

  if ((ptr->lat==NULL) || (ptr->lon==NULL) || 
      (ptr->mag==NULL) || (ptr->vertex==NULL) ||
      (ptr->ex==NULL) ||  (ptr->ey==NULL)) return -1;

  for (j=0;j<ngrid;j++) {
      ptr->vertex[4*poly]=num;
      ptr->vertex[4*num+1]=num-nlon+1;
      ptr->vertex[4*poly+2]=ptr->vertex[4*poly+1]+nlon;
      ptr->vertex[4*poly+3]=num+nlon;     
      poly++;  

      lat=(gptr[j].lat[0]+gptr[j].lat[1]+gptr[j].lat[2]+gptr[j].lat[3])/4;
      lon=(gptr[j].lon[0]+gptr[j].lon[1]+gptr[j].lon[2]+gptr[j].lon[3])/4;
      AACGMConvert(lat,lon,300.0,&mlat,&mlon,&r,0);

      ptr->lat[j]=mlat;
      ptr->lon[j]=mlon;
      ptr->mag[j]=0;
      num++;
  }

  ptr->num=ngrid;
  ptr->nlat=0;
  ptr->nlon=0;
  ptr->poly=poly;
  return 0;
}
 






