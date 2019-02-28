
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include "rtypes.h"
#include "rtime.h"
#include "rfile.h"
#include "aacgm.h"

#include "errstr.h"
#include "hlpstr.h"

#include "rmath.h"
#include "griddata.h"

#include "cnvgrid.h"
#include "cnvmap.h"
#include "cnvmapindex.h"
#include "cnvmapseek.h"
#include "cnvmapread.h"
#include "cnvmapsolve.h"
#include "make_pgrid.h"
#include "local_df_vel.h"


void shf_fit_vel(FILE *grdfp, time_t stime, struct ShfVel shf_vel){

  int yr,mo,dy,hr,mt;
  double sc;

  int i;

  int status;
 
  unsigned char help=0;
  unsigned char option=0;
 
  float latmin=50.0;
  float mlt;

  double glat,glon,r;

  int step=1;

  int num;
  int lat,lon;

  int *count=NULL;
   
  TimeEpochToYMDHMS(stime,&yr,&mo,&dy,&hr,&mt,&sc);
  rewind(grdfp);
  status=CnvMapFseek(grdfp,yr,mo,dy,hr,mt,sc,NULL,NULL);
  if (status ==-1) {
    fprintf(stderr,"File does not contain the requested interval.\n");
    exit(-1);
  }
  CnvMapFread(grdfp,map,grd);
 
  cgrid->type=1;

  CnvMapSolve(map,cgrid);

  for (i=0;i<cgrid->num;i++) {
    shf_vel[i].vx=cgrid->mag[i]*sind(grid->azm[i]);
    shf_vel[i].vy=cgrid->mag[i]*cosd(grid->azm[i]);
  }


}

