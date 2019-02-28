#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <zlib.h>
#include "rtypes.h"
#include "dmap.h"
#include "option.h"
#include "rtime.h"
#include "radar.h"
#include "rprm.h"
#include "rpos.h"
#include "fitdata.h"
#include "cfitdata.h"
#include "scandata.h"
#include "fitread.h"
#include "fitscan.h"
#include "fitindex.h"
#include "fitseek.h"
#include "rtypes.h"
#include "dmap.h"
#include "local_df_vel.h"
#include "invmag.h"
#include "griddata.h"

#include "cnvgrid.h"
#include "cnvmap.h"
#include "cnvmapindex.h"
#include "cnvmapseek.h"
#include "cnvmapread.h"
#include "cnvmapsolve.h"
#include "aacgmlib_v2.h"
#include "aacgm.h"
#include "igrflib.h"

#include <gsl/gsl_math.h> 
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define LENGTH(x,y) sqrt(x*x+y*y)

#define C 299792458.0 
#define PI 3.14159265359
#define LRE 6371.E3
#define MIN_RANGE 700
#define MAX(q,p) (((q)>(p))?(q):(p))
#define sind(x) (sin(fmod((x),360)*PI/180))
#define cosd(x) (cos(fmod((x),360)*PI/180))
#define tand(x) (tan(fmod((x),360)*PI/180))
#define MIN_ERR 1.
#define MAX_V 3500.
#define MAX_V_ERR 100.
#define MIN_V 30.
#define MIN_COUNT 3
#define MAX_BEAMS 24

float dlat;
float dlon;
long start_time;
long end_time;
int avg_ival;
int smooth;
double model_scale;
char *radar_list[30];
int nrad;
int nbp;
B_POINT* bp;
int ngrid;
int nedge;
CELL* grid;
NEIGHBOR* neighbors;
int *edge; 
M_ARRAY* edge_array;
double* edge_data;
M_ARRAY* kazm_array;
double* los_data;
double* los_err;
M_ARRAY* div_array;
M_ARRAY* div_array1;
M_ARRAY* smooth_array;
int mdays[]={31,59,90,120,151,181,212,243,273,304,334};
FILE *grd_file;
char shf_file[128];
void sub_sphazm();
int sub_sphcal();

int dayofweek(int d, int m, int y)
{
    static int t[] = { 0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4 };
    y -= m < 3;
    return ( y + y/4 - y/100 + y/400 + t[m-1] + d) % 7;
}
 
char *choppy( char *s )
{
  char *n = malloc( strlen( s ? s : "\n" ) );
  if( s )
    strcpy( n, s );
  if( s[strlen(s)-1] =='\n')
    n[strlen(n)-1]='\0';
  return n;
}


void parse_instructions(FILE *fp)
{
  char *line=NULL;
  size_t len=0;
  char *param=NULL;
  char *token;

  nrad=0;
  nbp=0;
  while( getline(&line,&len,fp)!= EOF ){
    if(line[0]=='#' || line[0]==' ')continue;
    param=strtok(line," ");
    if(strcmp(param,"boundary_point")==0)
      {
	nbp++;
	bp=realloc(bp,nbp*sizeof(B_POINT));
	sscanf(strtok(NULL," "),"%lf",&bp[nbp-1].lat);
	sscanf(strtok(NULL," "),"%lf",&bp[nbp-1].lon);
      }
   else if(strcmp(param,"lat_dl")==0)
     {
	sscanf(strtok(NULL," "),"%f",&dlat);
     }
   else if(strcmp(param,"lon_dl")==0)
     {
	sscanf(strtok(NULL," "),"%f",&dlon);
     }
   else if(strcmp(param,"start_time")==0)
     {
	sscanf(strtok(NULL," "),"%ld",&start_time);
     }
   else if(strcmp(param,"end_time")==0)
     {
	sscanf(strtok(NULL," "),"%ld",&end_time);
     }
   else if(strcmp(param,"avg_interval")==0)
     {
	sscanf(strtok(NULL," "),"%d",&avg_ival);
     }
   else if(strcmp(param,"shf_file")==0)
     {
	sscanf(strtok(NULL," "),"%s",shf_file);
     }
   else if(strcmp(param,"smooth")==0)
     {
	sscanf(strtok(NULL," "),"%d",&smooth);
     }
   else if(strcmp(param,"model_scale")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&model_scale);
     }
   else if(strcmp(param,"radar_list")==0)
     {
       while( (token=strtok(NULL," ")) != NULL)
	 {
	   radar_list[nrad]=choppy(token);
	   nrad++;
	 }
     }
 }
  if(line)
    free(line);
}

/* determines if a point is inside a polygon returns "1" for inside "0" for outside*/
int poly_in_out(double lat, double lon )
{
  int j,c=0;
  double lat1,lon1;
  double lat2,lon2;
  lat1=bp[0].lat;
  lon1=bp[0].lon;
  for( j=1; j<=nbp; j++)
    {
      lat2=bp[j % nbp].lat;
      lon2=bp[j % nbp].lon;
      if( fabs(lon2-lon1) > 180 )
	{
	  if( lon1 > 180 ) lon1-=360.;
	  if( lon2 > 180 ) lon2-=360.;
	}
      if( fabs(lon2-lon) > 180 ) lon-=360.;
      if( (((lat>=lat1) && (lat<=lat2)) || ((lat>=lat2) && (lat<=lat1))) &&
	  ( lon <(lon2-lon1)*(lat-lat1)/(lat2-lat1)+lon1) )
	c=!c;
      lat1=lat2;
      lon1=lon2;
    }
  return c;
}

/* determines if a point is inside a polygon returns "1" for inside "0" for outside*/
int cell_in_out(double lat, double lon, CELL cell )
{
  int j,c=0;
  double lat1,lon1;
  double lat2,lon2;
  lat1=cell.lat[0];
  lon1=cell.lon[0];
  for( j=1; j<=4; j++)
    {
      lat2=cell.lat[j % 4];
      lon2=cell.lon[j % 4];
      if( fabs(lon2-lon1) > 180 )
	{
	  if( lon1 > 180 ) lon1-=360.;
	  if( lon2 > 180 ) lon2-=360.;
	}
      if( fabs(lon2-lon) > 180 ) lon-=360.;
      if( (((lat>=lat1) && (lat<=lat2)) || ((lat>=lat2) && (lat<=lat1))) &&
	  ( lon <(lon2-lon1)*(lat-lat1)/(lat2-lat1)+lon1) )
	c=!c;
      lat1=lat2;
      lon1=lon2;
    }
  return c;
}

void make_box(double lat, double lon, double dx, double dy, CELL *box)
{
  double ldlat=180*dy/(LRE*PI);
  double dlon_low=180*dx/(LRE*cosd(lat-ldlat/2))/PI;
  double dlon_high=180*dx/(LRE*cosd(lat+ldlat/2))/PI;
  box->lat[0]=lat-ldlat/2;
  box->lat[1]=lat-ldlat/2;
  box->lat[2]=lat+ldlat/2;
  box->lat[3]=lat+ldlat/2;
  box->lon[0]=lon-dlon_low/2;
  box->lon[1]=lon+dlon_low/2;
  box->lon[2]=lon+dlon_high/2;
  box->lon[3]=lon-dlon_high/2;
}  


void set_boundaries(){
  
  double azm;
  double range;
  double lat,lon;
  double min_lat=91;
  double max_lat=-91;
  double min_lon=1000;
  double max_lon=-1000;
  void sub_gcp();
  int j,jp;

  for( jp=0; jp<nbp-1; jp++ )
    {
      if( bp[jp].lat < min_lat )min_lat=bp[jp].lat;
      if( bp[jp].lat > max_lat )max_lat=bp[jp].lat;
      if( bp[jp].lon < min_lon )min_lon=bp[jp].lon;
      if( bp[jp].lon > max_lon )max_lon=bp[jp].lon;
    }
  fprintf(grd_file,"min_lat: %f\n",min_lat); 
  fprintf(grd_file,"max_lat: %f\n",max_lat); 
  fprintf(grd_file,"min_lon: %f\n",min_lon);
  fprintf(grd_file,"max_lon: %f\n",max_lon);

  for( jp=0; jp<nbp-1; jp++ )
    {
      sub_sphazm(bp[jp].lon,bp[jp].lat,bp[jp+1].lon,bp[jp+1].lat,&azm,&range);
      double step=10;
      int nstep=range/step;
      for( j=1; j<=nstep; j++)
	{
	  range=(double)j*step;
	  sub_gcp(bp[jp].lon,bp[jp].lat,bp[jp+1].lon,bp[jp+1].lat,range,&lon,&lat);
	  fprintf(grd_file,"%f  %f\n",lat,lon);
	}
    }
  sub_sphazm(bp[nbp-1].lon,bp[nbp-1].lat,bp[0].lon,bp[0].lat,&azm,&range);
  double step=10;
  int nstep=range/step;
  for( j=1; j<=nstep; j++)
    {
      range=(double)j*step;
      sub_gcp(bp[nbp-1].lon,bp[nbp-1].lat,bp[0].lon,bp[0].lat,range,&lon,&lat);
      fprintf(grd_file,"%f  %f\n",lat,lon);
      if( lat < min_lat )min_lat=lat;
      if( lat > max_lat )max_lat=lat;
      if( lon < min_lon )min_lon=lon;
      if( lon > max_lon )max_lon=lon;
    }
  fprintf(grd_file,"0.0 0.0\n");
}

int edge_check(double lat, double lon, double dlon_l){
  int i_lat_m=poly_in_out(lat-dlat,lon);
  int i_lat_p=poly_in_out(lat+dlat,lon);
  int i_lon_m=poly_in_out(lat,lon+dlon_l);
  int i_lon_p=poly_in_out(lat,lon-dlon_l);
  return(!(i_lat_m & i_lat_p & i_lon_p & i_lon_m)); 
}

void make_grid()
{
  int j,jp;
  double lat,lon;
  double min_lat=91;
  double max_lat=-91;
  double min_lon=1000;
  double max_lon=-1000;

  for( jp=0; jp<nbp-1; jp++ )
    {
      if( bp[jp].lat < min_lat )min_lat=bp[jp].lat;
      if( bp[jp].lat > max_lat )max_lat=bp[jp].lat;
      if( bp[jp].lon < min_lon )min_lon=bp[jp].lon;
      if( bp[jp].lon > max_lon )max_lon=bp[jp].lon;
    }

  double start_lon,dlon_l;
  double center_lon=(max_lon+min_lon)/2;
  double center_lat=(max_lat+min_lat)/2;
  int nlon;
  int nlat=(max_lat-min_lat)/dlat+2;
  int jlat,jlon,eflg;

  double dx=fabs(dlon*DTOR*LRE*cosd(center_lat));
  double dy=fabs(dlat*DTOR*LRE);

  /*    
    walks through a latitude-longitude grid and determines interior points
    prints the latitude and longitude of the interior points
    prints lats & lons of corners of dx x dy box surrounding ponits
   */

  ngrid=0;
  nedge=0;
  fprintf(stderr,"MAKE GRID: min_lat=%f max_lat=%f\n",min_lat,max_lat);
  for( jlat=0; jlat<nlat; jlat++ )
    {
      lat=min_lat+((double)jlat+.5)*dlat;
      dlon_l=180*dx/(LRE*cosd(lat))/PI;
      nlon=(int)(max_lon-min_lon)/dlon_l;
      start_lon=center_lon-(double)nlon*dlon_l/2.;
      for( jlon=0; jlon<nlon; jlon++)
	{
	  lon=start_lon+((double)jlon+.5)*dlon_l;
	  if(poly_in_out(lat,lon)) 
	    {
	      if(edge_check(lat,lon,dlon_l)){
		edge=realloc(edge,(ngrid+1)*sizeof(int));
		edge[nedge]=ngrid;
		nedge++;
	      }
	      fprintf(grd_file,"%f  %f %d\n",lat,lon,edge_check(lat,lon,dlon_l));
	      grid=realloc(grid,(ngrid+1)*sizeof(CELL));
	      make_box(lat,lon,dx,dy,&grid[ngrid]);
	      fprintf(grd_file,"%f  %f  %f  %f\n",grid[ngrid].lat[0],grid[ngrid].lat[1],grid[ngrid].lat[2],grid[ngrid].lat[3]);
	      fprintf(grd_file,"%f  %f  %f  %f\n",grid[ngrid].lon[0],grid[ngrid].lon[1],grid[ngrid].lon[2],grid[ngrid].lon[3]);
	      ngrid++;
	    }
	}
    }  
  fprintf(grd_file,"0.0 0.0 0.0 \n");
  fprintf(stderr,"NGRID = %d\n",ngrid);
}

void find_neighbors(){
  int jg,ig;
  int jc;
  double lat,lon;
  for( jg=0; jg<ngrid; jg++ ){
    for( jc=0; jc<4; jc++){
      if(jc <= 1) lat=grid[jg].lat[jc]-.5*dlat; else lat=grid[jg].lat[jc]+.5*dlat;
      lon=grid[jg].lon[jc];
      neighbors[jg].c[jc]=-1;
      for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig]))neighbors[jg].c[jc]=ig;
    }
  }
}

void set_edge_ar(){
  int je;
  double lat,lon;
  int j,jp;
  double min_lat=91;
  double max_lat=-91;
  double min_lon=1000;
  double max_lon=-1000;
  for( jp=0; jp<nbp-1; jp++ )
    {
      if( bp[jp].lat < min_lat )min_lat=bp[jp].lat;
      if( bp[jp].lat > max_lat )max_lat=bp[jp].lat;
      if( bp[jp].lon < min_lon )min_lon=bp[jp].lon;
      if( bp[jp].lon > max_lon )max_lon=bp[jp].lon;
    }
  double* coef;
  double center_lon=(max_lon+min_lon)/2;
  double center_lat=(max_lat+min_lat)/2;

  edge_array=calloc(2*nedge,sizeof(struct mod_array));

  for( je=0; je<nedge; je++){
    edge_array[2*je].coef=calloc(2*ngrid,sizeof(double));
    edge_array[2*je+1].coef=calloc(2*ngrid,sizeof(double));
    edge_array[2*je].coef[edge[je]]=1;
    edge_array[2*je+1].coef[edge[je]+ngrid]=1;
  }
}

void make_div_ar(){
  int igr,iedg,inr0,inr1;
  double fr0,fr1;
  int idv=0;
  int flg;
  double dx,dlon_l;
  double dy=fabs(dlat*DTOR*LRE);

  div_array=calloc(ngrid-nedge,sizeof(struct mod_array));
  for( igr=0; igr<ngrid; igr++){
    flg=0;
    for(iedg=0; iedg<nedge; iedg++) if( edge[iedg]==igr )flg=1;
    if( flg==0 ){
      div_array[idv].coef=calloc(2*ngrid,sizeof(double));
      dlon_l=fabs(grid[igr].lon[1]-grid[igr].lon[0]);
      dx=fabs(dlon_l*DTOR*LRE*cosd((grid[igr].lat[0]+grid[igr].lat[2])/2));
      div_array[idv].coef[igr-1]=-1/dx;
      div_array[idv].coef[igr]=1/dx;
      /*      determine neighbor lengths */
      div_array[idv].coef[igr+ngrid]=1/dy;
      inr0=neighbors[igr].c[0];
      inr1=neighbors[igr].c[1];
      if( inr0 == -1 ){ 
	div_array[idv].coef[inr1+ngrid]=-1/dy;
      }else if(inr1 == -1){
	div_array[idv].coef[inr0+ngrid]=-1/dy;
      } else { 
	fr0=fabs((grid[inr0].lon[2]-grid[igr].lon[0])/(grid[igr].lon[1]-grid[igr].lon[0]));
	fr1=1-fr0;
	div_array[idv].coef[inr0+ngrid]=-fr0/dy;
	div_array[idv].coef[inr1+ngrid]=-fr1/dy;
      }
      idv++;
    }
  }
}

void make_div_ar1(){
  int igr,iedg,inr0,inr1;
  double fr0,fr1;
  int idv=0;
  int flg;
  double dx,dlon_l;
  double dy=fabs(dlat*DTOR*LRE);

  div_array1=calloc(ngrid-nedge,sizeof(struct mod_array));
  for( igr=0; igr<ngrid; igr++){
    flg=0;
    for(iedg=0; iedg<nedge; iedg++) if( edge[iedg]==igr )flg=1;
    if( flg==0 ){
      div_array1[idv].coef=calloc(2*ngrid,sizeof(double));
      dlon_l=fabs(grid[igr].lon[1]-grid[igr].lon[0]);
      dx=fabs(dlon_l*DTOR*LRE*cosd((grid[igr].lat[0]+grid[igr].lat[2])/2));
      div_array1[idv].coef[igr+1]=1/dx;
      div_array1[idv].coef[igr]=-1/dx;
      /*      determine neighbor lengths */
      div_array1[idv].coef[igr+ngrid]=-1/dy;
      inr0=neighbors[igr].c[2];
      inr1=neighbors[igr].c[3];
      if( inr0 == -1 ){ 
	div_array1[idv].coef[inr1+ngrid]=1/dy;
      }else if(inr1 == -1){
	div_array1[idv].coef[inr0+ngrid]=1/dy;
      } else { 
	fr0=fabs((grid[inr0].lon[2]-grid[igr].lon[0])/(grid[igr].lon[1]-grid[igr].lon[0]));
	fr1=1-fr0;
	div_array1[idv].coef[inr0+ngrid]=fr0/dy;
	div_array1[idv].coef[inr1+ngrid]=fr1/dy;
      }
      idv++;
    }
  }
}

void make_smooth_ar(){
  int igr,iedg,inr0,inr1;
  double fr0,fr1;
  int ism=0;
  double count;
  int flg;
  double dx,dlon_l;
  double dy=fabs(dlat*DTOR*LRE);

  smooth_array=calloc(2*(ngrid-nedge),sizeof(struct mod_array));
  for( igr=0; igr<ngrid; igr++){
    flg=0;
    count=0;
    for(iedg=0; iedg<nedge; iedg++) if( edge[iedg]==igr )flg=1;
    if( flg==0 ){
      smooth_array[ism].coef=calloc(2*ngrid,sizeof(double));
      smooth_array[ism+1].coef=calloc(2*ngrid,sizeof(double));
       
      inr0=neighbors[igr].c[0];
      if( inr0 != -1){
	smooth_array[ism].coef[inr0]=-1.;       /* vx for south west neighbor */
	smooth_array[ism+1].coef[inr0+ngrid]=-1.; /* vy for south west neighbor */
	count+=1.;
      }	
      inr0=neighbors[igr].c[1];
      if( inr0 != -1){
	smooth_array[ism].coef[inr0]=-1.;       /* vx for south east neighbor */
	smooth_array[ism+1].coef[inr0+ngrid]=-1.; /* vy for south east neighbor */
	count+=1.;
      }	      
      inr0=neighbors[igr].c[2];
      if( inr0 != -1){
	smooth_array[ism].coef[inr0]=-1.;       /* vx for north east neighbor */
	smooth_array[ism+1].coef[inr0+ngrid]=-1.; /* vy for north east neighbor */
	count+=1.;
      }	
      inr0=neighbors[igr].c[3];
      if( inr0 != -1){
	smooth_array[ism].coef[inr0]=-1.;       /* vx for south east neighbor */
	smooth_array[ism+1].coef[inr0+ngrid]=-1.; /* vy for south east neighbor */
	count+=1.;

      }	
      smooth_array[ism].coef[igr-1]=-1; /* vx for west neighbor */
      smooth_array[ism+1].coef[igr-1+ngrid]=-1; /* vy for west neighbor */
      count+=1.;
      smooth_array[ism].coef[igr+1]=-1; /* vx for west neighbor */
      smooth_array[ism+1].coef[igr+1+ngrid]=-1; /* vy for west neighbor */
      count+=1.;
      for( inr0=0; inr0<2*ngrid; inr0++)smooth_array[ism].coef[inr0]/=count;
      for( inr0=0; inr0<2*ngrid; inr0++)smooth_array[ism+1].coef[inr0]/=count;

      smooth_array[ism].coef[igr]=1; /* vx for grid point */      
      smooth_array[ism+1].coef[igr+ngrid]=1; /* vy for grid point */
      ism+=2;
    }
  }

}



struct tm *parse_date_str( long t_i)
{
  struct tm *t_o;
  char *tz;

  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();

  t_o=malloc(sizeof(struct tm));
  t_o->tm_year=(int)(t_i/1e8);
  t_o->tm_mon=(int)((t_i-1e8*t_o->tm_year)/1e6);
  t_o->tm_mday=(int)((t_i-1e8*t_o->tm_year-1e6*t_o->tm_mon)/1e4);
  t_o->tm_hour=(int)((t_i-1e8*t_o->tm_year-1e6*t_o->tm_mon-t_o->tm_mday*1e4)/1e2);
  t_o->tm_min=(int)(t_i-1e8*t_o->tm_year-1e6*t_o->tm_mon-t_o->tm_mday*1e4-t_o->tm_hour*1e2);
  t_o->tm_yday=mdays[t_o->tm_mon-1]+t_o->tm_mday;
  t_o->tm_wday=dayofweek(t_o->tm_mday,t_o->tm_mon,t_o->tm_year);
  if( IS_LEAPYEAR(t_o->tm_year) && t_o->tm_mon>2) t_o->tm_yday++;
  t_o->tm_sec=0;
  t_o->tm_year-=1900;
  t_o->tm_mon-=1;
  t_o->tm_isdst=0;
  return t_o;
}

time_t fname_to_time(char *fname)
{
  int datev;
  char datestr[10]="0";
  char hr_str[3]="0";
  char mn_str[3]="0";
  struct tm f_tm;
  char *tz;
  
  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();
  
  memset(datestr, '\0', sizeof datestr);
  memset(hr_str, '\0', sizeof hr_str);
  memset(mn_str, '\0', sizeof mn_str);

  strncpy(datestr,fname,8);
  strncpy(hr_str,fname+9,2);
  strncpy(mn_str,fname+11,2);
  datev=atoi(datestr);
  f_tm.tm_year=datev/10000;
  f_tm.tm_mon=(datev-10000*f_tm.tm_year)/100;
  f_tm.tm_mday=datev-10000*f_tm.tm_year-100*f_tm.tm_mon;
  f_tm.tm_yday=mdays[f_tm.tm_mon-1]+f_tm.tm_mday;
  if( IS_LEAPYEAR(f_tm.tm_year) && f_tm.tm_mon>2) f_tm.tm_yday++;
  f_tm.tm_hour=atoi(hr_str);
  f_tm.tm_min=atoi(mn_str);
  f_tm.tm_sec=0;
  f_tm.tm_isdst=0;
  f_tm.tm_wday=dayofweek(f_tm.tm_mday,f_tm.tm_mon,f_tm.tm_year);
  f_tm.tm_year-=1900; /* unix epoch year correction */
  f_tm.tm_mon-=1;     /* unix epoch month 0 to 11 */

  return mktime(&f_tm);
}

FILE_INFO *file=NULL;
FILE_INFO *select_file(char *radar, time_t time)
{
  char yr_str[5], mo_str[3], dy_str[3];
  char *raid_path;
  char dir_path[PATH_LEN];
  time_t ftime;
  time_t diftime;
  time_t mindif=100000;
  DIR *dp;
  struct dirent *ep;
  char *tz;

  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();

  raid_path=getenv("RAID_PATH");
  struct tm *in_time;
  in_time=gmtime(&time);
  int yr=in_time->tm_year+1900;
  int mo=in_time->tm_mon+1;
  int dy=in_time->tm_mday;
  int hr=in_time->tm_hour;
  int mn=in_time->tm_min;
  int sc=in_time->tm_sec;

  CNV_TO_STR(yr,yr_str);
  CNV_TO_STR(mo,mo_str);
  CNV_TO_STR(dy,dy_str);
  
  sprintf(file->dir_path,"%s","");
  sprintf(file->fname,"%s","");

  sprintf(dir_path,"%s%s/%s.%s",raid_path,yr_str,mo_str,dy_str);
  fprintf(stderr,"\n %s%s/%s.%s\n",raid_path,yr_str,mo_str,dy_str);
  if((dp=opendir(dir_path))==NULL) {
    fprintf(stderr,"---- COULDN'T OPEN DATA DIRECTORY ----\n");
    return(file);
  }

  while( (ep=readdir(dp)) )
    {
      if( strstr(ep->d_name,".gz") != NULL) continue;
      if( strstr(ep->d_name,".bz2") != NULL) continue;
      if( strstr(ep->d_name,radar) != NULL)
	{
	  ftime=fname_to_time(ep->d_name);
	  diftime=time-ftime;
	  if( diftime>=0 && diftime<mindif )
	  {
	    sprintf(file->dir_path,"%s",dir_path);
	    sprintf(file->fname,"%s",ep->d_name);
	    mindif=diftime;
	  }
	}
    }
  return file;
}

time_t ftime(struct RadarParm *prm){
  int yr,mo,dy,hr,mt,sc;
  yr=prm->time.yr;
  mo=prm->time.mo;
  dy=prm->time.dy;
  hr=prm->time.hr;
  mt=prm->time.mt;
  sc=prm->time.sc;
  return (time_t)TimeYMDHMSToEpoch(yr,mo,dy,hr,mt,sc);
}


void get_pos_ar(struct RadarSite *site, struct RadarParm *prm, struct RadarPos *rdrpos){
  
  double rho,lat,lon;
  double geoazm,elv;
  int rn,bm,rsep,frang;
  int rxrise,yr;
  int chis=1;
  rsep=prm->rsep;
  frang=prm->frang;
  rxrise=prm->rxrise;
  yr=prm->time.yr;
  
  for( bm=0; bm<site->maxbeam; bm++) for( rn=0; rn<=prm->nrang; rn++){
      
      RPosGeo(0,bm,rn,site,frang,rsep,
	      site->recrise,0,&rho,&lat,&lon,chis);
      rdrpos->lat[bm][rn]=lat;
      rdrpos->lon[bm][rn]=lon;
      RPosRngBmAzmElv(bm,rn,yr,site,frang,rsep,rxrise,300.0,&geoazm,&elv,chis);
      rdrpos->kazm[bm][rn]=geoazm;
    }
}

void map_pos_to_grid(struct RadarPos rdrpos, struct RadarMap *rmap){
  int ib,ir,ig;
  double lat,lon;
  for( ib=0; ib<MAX_BEAM; ib++)for( ir=0; ir<MAX_RANGE; ir++){
      lat=rdrpos.lat[ib][ir];
      lon=rdrpos.lon[ib][ir];
      rmap->cell[ib][ir]=-1;
      for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig]))rmap->cell[ib][ir]=ig;
    }
}

double grid_lat(int ig){
  return((grid[ig].lat[0]+grid[ig].lat[1]+grid[ig].lat[2]+grid[ig].lat[3])/4);
}

double grid_lon(int ig){
  return((grid[ig].lon[0]+grid[ig].lon[1]+grid[ig].lon[2]+grid[ig].lon[3])/4);
}


void convert_vec_to_geo(double lat, double lon, double mag, double azm_in, double *azm_out){

  double lat1,lon1;
  double glat,glon,r;
  double glat1,glon1;
  int stat;

  stat=sub_sphcal(lon,lat,azm_in,mag,&lon1,&lat1);
  
  stat=AACGM_v2_Convert(lat,lon,300, &glat, &glon, &r,1);
  stat=AACGM_v2_Convert(lat1,lon1,300, &glat1, &glon1, &r,1);

  sub_sphazm(glon,glat,glon1,glat1,azm_out,&r);


}

struct CnvMapData *map;
struct GridData *grd;
struct CnvGrid *cgrid;

void shf_fit_vel(FILE *grdfp, time_t stime, struct ShfVel *shf_vel, double *shf_error){

  int yr,mo,dy,hr,mt;
  double sc;
  int i;
  int old=0;
  int status;
  double az_out;
  FILE *shfvf;
  
  TimeEpochToYMDHMS(stime,&yr,&mo,&dy,&hr,&mt,&sc);
  fprintf(stderr,"SHF; %d : %d : %d\n",hr,mt,(int) sc);
  rewind(grdfp);
  status=CnvMapFseek(grdfp,yr,mo,dy,hr,mt,sc,NULL,NULL);
  if (status ==-1) {
    fprintf(stderr,"File does not contain the requested interval.\n");
    exit(-1);
  }
  CnvMapFread(grdfp,map,grd);

  IGRF_SetDateTime(yr,mo,dy,hr,mt,(int)sc);
  AACGM_v2_SetDateTime(yr,mo,dy,hr,mt,(int)sc);
  
  float decyear=(float)yr+(float)mo/12+(float)dy/30;   

  cgrid->type=1;
  CnvMapSolve(map,cgrid,decyear,old);
  fprintf(stderr,"shf error %f\n",map->rms_err);
  *shf_error=map->rms_err;
  /* shfvf=fopen("shfvel_out","w"); */
  /* fprintf(shfvf,"%d\n",cgrid->num); */
  for (i=0;i<cgrid->num;i++) {
    
    convert_vec_to_geo(cgrid->lat[i],cgrid->lon[i],cgrid->mag[i],cgrid->azm[i],&az_out);
    if( *shf_error >= 1.e3 )*shf_error=1.e3;

    shf_vel[i].vx=cgrid->mag[i]*sin(az_out*(double)DTOR);
    shf_vel[i].vy=cgrid->mag[i]*cos(az_out*(double)DTOR);
    
     /* fprintf(shfvf,"%f %f %f %f %f\n",cgrid->lat[i],cgrid->lon[i],cgrid->mag[i],cgrid->mag[i],az_out); */
        /* fprintf(shfvf,"%f %f %f %f %f\n",cgrid->lat[i],cgrid->lon[i],shf_vel[i].vx,shf_vel[i].vy,az_out); */
    
      /*    
    if( *shf_error < 1.e5 ){ 
      shf_vel[i].vx=cgrid->mag[i]*sind(az_out);
      shf_vel[i].vy=cgrid->mag[i]*cosd(az_out);
    }else{
      *shf_error=1.e5;
      shf_vel[i].vx=50.;
      shf_vel[i].vy=50.;
    }
    */
  }
  /* fclose(shfvf); */

}


struct RadarParm *prm;
struct FitData *fit;
struct FitIndex *inx;

struct RadarNetwork *network;  
struct Radar *radar;

struct RadarPos *rpos;

extern int pinv(int m,int n,double** a,double** inv);

	       
int main(int argc, char *argv[]){

  /* read instruction file
     file should have date, time, lat and lon of area corners, radars to contribute
     filtering instructions
  */
  FILE *fp;
  FILE *vf;
  FILE *shfvf;
  FILE *rptrs[NRAD];
  FILE *grdfp;
  char file_name[128];
  char *flist[MXFILES];
  int jj,nf,colcount;
  int nbms;
  struct tm *t_start;
  struct tm *t_end; 
  time_t time, ssec, esec, istart, iend;
  char *envstr;
  struct RadarSite *site=NULL;
  struct RadarMap *rmap;
  struct ShfVel *shf_vel;
  int *grid_data_count;
  gsl_matrix *coef;
  gsl_matrix *coefTrans;
  gsl_vector *data;
  gsl_vector *GM;
  gsl_matrix *CD;  
  gsl_vector *solution;
  gsl_matrix *result;

  gsl_matrix * inv;

  double **aa;
  double **pinv_ar;


  smooth=0;
  model_scale=100.;

  if (argc <= 1){
    fprintf(stderr,"********NO INSTRUCTION FILE GIVEN********\n");
    exit(-1);
  }
  strcpy(file_name,argv[1]);
  if ((fp=fopen(file_name,"r")) == 0) {
    fprintf(stderr,"******FIT INSTRUCTION FILE %s NOT FOUND******\n",argv[1]);
    exit(-1);
  }
  /* parse file */
  parse_instructions(fp);
  fclose(fp);

  /* create local grid
   */
  grd_file=fopen("grid.dat","w");
  set_boundaries();
  make_grid();
  fclose(grd_file);

  neighbors=calloc(ngrid,sizeof(struct Neighbor));
  find_neighbors();

  /* set boundary vel over array */
  set_edge_ar();

  /* set divergence = 0 over array */
  make_div_ar();
  make_div_ar1();

  /* set smoothing prior */
  if(smooth)make_smooth_ar();

  envstr=getenv("SD_RADAR");
  if (envstr==NULL) {
    fprintf(stderr,"Environment variable 'SD_RADAR' must be defined.\n");
    exit(-1);
  }

  fp=fopen(envstr,"r");

  if (fp==NULL) {
    fprintf(stderr,"Could not locate radar information file.\n");
    exit(-1);
  }

  network=RadarLoad(fp);
  fclose(fp); 
  if (network==NULL) {
    fprintf(stderr,"Failed to read radar information.\n");
    exit(-1);
  }

  envstr=getenv("SD_HDWPATH");
  if (envstr==NULL) {
    fprintf(stderr,"Environment variable 'SD_HDWPATH' must be defined.\n");
    exit(-1);
  }

  RadarLoadHardware(envstr,network);

  /* find fit files and read into model-array and data vector
     requires determining k-vectors
   */
  t_start=parse_date_str(start_time);
  t_end=parse_date_str(end_time);

  time=mktime(t_start);
  esec=mktime(t_end);  


  if(file)free(file);
  file=calloc(1,sizeof(FILE_INFO));

  nrad--;
  for( jj=0; jj<nrad; jj++ )
    {
      rptrs[jj]=NULL;
      file=select_file(radar_list[jj],time);
      if ( file->fname !=NULL){
	fprintf(stderr,"%s%s\n",file->dir_path,file->fname);
	sprintf(file_name,"%s/%s",file->dir_path,file->fname);
	if(rptrs[jj])fclose(rptrs[jj]);
	if((rptrs[jj]=fopen(file_name,"r"))== NULL)fprintf(stderr,"could not open file: %s",file_name);
      }
    }
  fprintf(stderr,"first files open\n");
  int yr,mo,dy,hr,mt;
  int yr_st,mo_st,dy_st,hr_st,mt_st;
  double sc,sc_st;
  int s,jr,jrmin,jc;
  int jjc,jjr,count;
  int jcmn,jcmx,jrmn,jrmx;
  int frang,rsep;
  int *rsepl;
  int *frangl;
  int *nrangl;
  double kazm;
  double avg_val;
  double dev_val;
  int ndata=0;
  int grid_cell;
  int center=1;
  time_t filetime;
  int neqn;
  int rval;
  int jeqn;

  double ModelSigma=1.e3;
  double shf_error;
  double divscale;
  double temp;

  double **filt_ar;
  double **filt_ar1;
  double **filt_dev;
  double **filt_err;
  double **filt_cnt;
  int **good_ar;
  
  int m, n, lda, ldu, ldvt, info, lwork;
  double wkopt;
  double* work;
  /* Local arrays */
  /*
  double s[N], u[LDU*M], vt[LDVT*N];
  double a[LDA*N]; 
*/
  /*  double *s, *u, *vt;*/
  double *a; 
  clock_t t1,t2;


  nrangl=calloc(nrad, sizeof(int));
  rsepl=calloc(nrad, sizeof(int));
  frangl=calloc(nrad, sizeof(int));

  prm=RadarParmMake();
  fit=FitMake();
  rpos=calloc(nrad, sizeof(struct RadarPos));
  rmap=calloc(nrad, sizeof(struct RadarMap));
  

  /* set dummy vel over interior of array */

  /* fill in data array */

  shf_vel=calloc(ngrid, sizeof(struct ShfVel));
  grid_data_count=malloc(ngrid*sizeof(int));

  grdfp=fopen(shf_file,"r");
  if (grdfp==NULL) {
     fprintf(stderr,"File not found.\n");
     exit(-1);
  }


  grd=GridMake();
  map=CnvMapMake();
  cgrid=CnvGridMake();
  make_local_pgrid(grid,ngrid,cgrid);

  if (CnvMapFread(grdfp,map,grd)==-1) {
    fprintf(stderr,"Error reading shf file.\n");
    exit(-1);
  }

  vf=fopen("vel_out","w");
  fp=fopen("coef_array","w");
  unsigned long position;


  while( time<esec )
    {	
      TimeEpochToYMDHMS(time-avg_ival/2,&yr_st,&mo_st,&dy_st,&hr_st,&mt_st,&sc_st);
      printf("%d %d %d %d %d %f\n",yr_st,mo_st,dy_st,hr_st,mt_st,sc_st);
      ndata=0;
      for(jj=0; jj<ngrid; jj++)grid_data_count[jj]=0;
      for( jj=0; jj<nrad; jj++){


	if( rptrs[jj]!=NULL ){
	  s=fseek(rptrs[jj],0L,SEEK_SET);
	  if( s==0 )s=FitFseek(rptrs[jj],yr_st,mo_st,dy_st,hr_st,mt_st,sc_st,NULL,inx);
	  if( s==0 )s=FitFread(rptrs[jj],prm,fit);
	  if( s==0 )fprintf(stderr,"%s file time: %d %d %d %d %d %d %d %ld\n",radar_list[jj],prm->time.yr,
			    prm->time.mo,prm->time.dy,prm->time.hr,prm->time.mt,prm->time.sc,prm->stid,
			    ftime(prm)-time);
	  if( s!=0 ){
	    fprintf(stderr,"File does not contain the requested interval. %d:%d\n",hr_st,mt_st);
	    if(rptrs[jj]){
	      fprintf(stderr,"*******Closing %s file\n",radar_list[jj]);
	      fclose(rptrs[jj]);
	      rptrs[jj]=NULL;
	    }
	  }
	}

	if( rptrs[jj]==NULL ){
	  file=select_file(radar_list[jj],time-avg_ival);
	  if( file->fname==NULL ) continue;
	  sprintf(file_name,"%s/%s",file->dir_path,file->fname);
	  fprintf(stderr,"*******Opening file %s\n",file_name); 
	  rptrs[jj]=fopen(file_name,"r");
	  if( rptrs[jj]!=NULL ){
	    s=FitFseek(rptrs[jj],yr_st,mo_st,dy_st,hr_st,mt_st,sc_st,NULL,inx);
	    if( s==0 )s=FitFread(rptrs[jj],prm,fit);
	    TimeEpochToYMDHMS(time+avg_ival,&yr,&mo,&dy,&hr,&mt,&sc);
	  }
	}
	  
	radar=RadarGetRadar(network,prm->stid);
	site=RadarYMDHMSGetSite(radar,yr_st,mo_st,dy_st,hr_st,mt_st,(int) 0);
	if (site==NULL) {fprintf(stderr,"NULL site\n"); continue;}

	nbms=(int)site->maxbeam;
	filt_ar=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_ar[jr]=(double*)calloc(prm->nrang,sizeof(double));
	filt_ar1=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_ar1[jr]=(double*)calloc(prm->nrang,sizeof(double));
	filt_dev=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_dev[jr]=(double*)calloc(prm->nrang,sizeof(double));
	filt_err=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_err[jr]=(double*)calloc(prm->nrang,sizeof(double));
	filt_cnt=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_cnt[jr]=(double*)calloc(prm->nrang,sizeof(double));
	good_ar=malloc(nbms*sizeof(int*));
	for( jr=0; jr<nbms; jr++)good_ar[jr]=(int*)calloc(prm->nrang,sizeof(int));
	nrangl[jj]=prm->nrang;

	if( rptrs[jj]!=NULL ){
	  while( ftime(prm)<time+avg_ival/2 ){
	    if( s!=0 ){
	      fprintf(stderr,"File does not contain the requested interval. %d:%d\n",hr_st,mt_st);
	      filetime=ftime(prm);
	      fprintf(stderr,"%s  %ld\n",radar_list[jj],time+avg_ival);
	      file=select_file(radar_list[jj],time+avg_ival);
	      sprintf(file_name,"%s/%s",file->dir_path,file->fname);
	      if(rptrs[jj]){
		fprintf(stderr,"*******Closing %s file\n",radar_list[jj]);
		fclose(rptrs[jj]);
	      }
	      rptrs[jj]=NULL;
	      if( file->fname == NULL )continue;
	      fprintf(stderr,"*******Opening file %s\n",file_name); 
	      if( (rptrs[jj]=fopen(file_name,"r")) ==NULL )continue;
	      /*	      rewind(rptrs[jj]);*/
	      TimeEpochToYMDHMS(time+avg_ival,&yr,&mo,&dy,&hr,&mt,&sc);
	      if((s=FitFseek(rptrs[jj],yr,mo,dy,hr,mt,sc,NULL,inx))!=0 ||
		 (s=FitFread(rptrs[jj],prm,fit))!=0){
		fprintf(stderr,"problem reading %s\n",file_name);
		fclose(rptrs[jj]);
		rptrs[jj]=NULL;
		break;
	      }
	      fprintf(stderr,"main: %s\n",file->fname);
	    }


	    if( s==0 ){

	      rsep=prm->rsep;
	      frang=prm->frang;
	      TimeEpochToYMDHMS(time,&yr,&mo,&dy,&hr,&mt,&sc);
	      if((rsep != rsepl[jj]) || (frang != frangl[jj])){
		fprintf(stderr,"recalculating position array %s %d %d  %d  %d\n",
			radar_list[jj],rsep,rsepl[jj],frang,frangl[jj]);
		radar=RadarGetRadar(network,prm->stid);
		site=RadarYMDHMSGetSite(radar,yr,mo,dy,hr,mt,(int) sc);
		if (site==NULL){fprintf(stderr,"site was null\n"); continue;}
		get_pos_ar(site,prm,&rpos[jj]);
		map_pos_to_grid(rpos[jj],&rmap[jj]);
		rsepl[jj]=rsep;
		frangl[jj]=frang;
	      }
	      if(prm->nrang != nrangl[jj]){
		fprintf(stderr,"%d station: %s  nrang mismatch old: %d   new: %d...reallocating\n",jj,radar_list[jj],nrangl[jj],prm->nrang);

		for( jr=0; jr<nbms; jr++)if(filt_ar[jr] != NULL)free(filt_ar[jr]); 
		for( jr=0; jr<nbms; jr++)if(filt_ar1[jr] != NULL)free(filt_ar1[jr]); 
		for( jr=0; jr<nbms; jr++)if(filt_dev[jr] != NULL)free(filt_dev[jr]); 
		for( jr=0; jr<nbms; jr++)if(filt_err[jr] != NULL)free(filt_err[jr]); 
		for( jr=0; jr<nbms; jr++)if(filt_cnt[jr] != NULL)free(filt_cnt[jr]); 
		for( jr=0; jr<nbms; jr++)if(good_ar[jr] != NULL)free(good_ar[jr]);


		for( jr=0; jr<nbms; jr++)filt_ar[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)filt_ar1[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)filt_dev[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)filt_err[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)filt_cnt[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)good_ar[jr]=(int*)calloc(prm->nrang,sizeof(int));
		
		if(prm->nrang == 0)for( jr=0; jr<nbms; jr++){
		    filt_ar[jr]=NULL;
		    filt_ar1[jr]=NULL;
		    filt_dev[jr]=NULL;
		    filt_err[jr]=NULL;
		    filt_cnt[jr]=NULL;
		    good_ar[jr]=NULL;
		  }		  
		nrangl[jj]=prm->nrang;
	      }
	      jrmin=(MIN_RANGE-frang)/rsep;
	      if( fabs(ftime(prm)-time) < (double)avg_ival/2){
		yr=prm->time.yr;
		mo=prm->time.mo;
		dy=prm->time.dy;
		hr=prm->time.hr;
		mt=prm->time.mt;
		sc=(double)prm->time.sc;
		for( jr=jrmin; jr<prm->nrang; jr++ )
		  if(fit->rng[jr].qflg == 1 && fit->rng[jr].gsct == 0 && 
		     (fit->rng[jr-1].qflg == 1 || fit->rng[jr+1].qflg == 1) &&
		     (fit->rng[jr-1].gsct == 0 || fit->rng[jr+1].gsct == 0) &&
		     fabs(fit->rng[jr].v_err) < MAX_V_ERR &&
		     fabs(fit->rng[jr].v) < MAX_V && fabs(fit->rng[jr].v) > MIN_V){
		    
		    filt_ar[prm->bmnum][jr]+=fit->rng[jr].v;
		    filt_err[prm->bmnum][jr]+=fit->rng[jr].v_err;
		    filt_cnt[prm->bmnum][jr]+=1.;
		    
		  }
	      }
	    }
	    s=FitFread(rptrs[jj],prm,fit);
	  }
	}

	for( jc=0; jc<site->maxbeam; jc++ )for( jr=1; jr<prm->nrang; jr++){
	    jcmn=MAX(jc-1,0);
	    jcmx=MIN(jc+1,site->maxbeam-1);
	    jrmn=MAX(jr-1,0);
	    jrmx=MIN(jr+1,prm->nrang-1);
	    count=0;
	    avg_val=0;
	    dev_val=0;
	    for(jjc=jcmn; jjc<=jcmx; jjc++)for(jjr=jrmn; jjr<=jrmx; jjr++){
		if( filt_ar[jjc][jjr]!=0 ){
		  count++;
		  avg_val+=(filt_ar[jjc][jjr]/filt_cnt[jjc][jjr]);
		}
	      }
	    //	    fprintf(stderr,"%d %d %d %8.2f",jc,jr,count,avg_val/(double)count) ;
	    if(count<MIN_COUNT){
	      good_ar[jc][jr]=0;
	      //	      fprintf(stderr,"\n");
	      continue;
	    }
	    avg_val/=(double)count;
	    for(jjc=jcmn; jjc<=jcmx; jjc++)for(jjr=jrmn; jjr<=jrmx; jjr++){
		if( filt_ar[jjc][jjr]!=0 )
		  dev_val+=(filt_ar[jjc][jjr]-avg_val)*(filt_ar[jjc][jjr]-avg_val)/(filt_cnt[jjc][jjr]*filt_cnt[jjc][jjr]);	
	      }
	    good_ar[jc][jr]=1;
	    filt_ar1[jc][jr]=avg_val;
	    filt_dev[jc][jr]=sqrt(dev_val/(double)count);
	    //	    fprintf(stderr,"%8.2f\n",filt_dev[jc][jr]);
	    //	    fprintf(stderr,"%d %d  %8.2f %8.2f %8.2f\n",jc,jr,filt_ar[jc][jr],filt_ar1[jc][jr],filt_dev[jc][jr]);
	  }
	for( jc=0; jc<nbms; jc++ )for( jr=1; jr<prm->nrang; jr++){
	    if ((grid_cell=rmap[jj].cell[jc][jr]) != -1 && good_ar[jc][jr] == 1 &&
		filt_ar[jc][jr] > filt_ar1[jc][jr]-3*filt_dev[jc][jr] &&
		filt_ar[jc][jr] < filt_ar1[jc][jr]+3*filt_dev[jc][jr] && filt_cnt[jc][jr]>0 ){		 		    

	      grid_data_count[grid_cell]++;
		    
	      kazm=rpos[jj].kazm[jc][jr];
		    
	      printf("%d %d %d %f %f %f %f %f %f %f\n",prm->stid,jc,jr,rpos[jj].lat[jc][jr],
		     rpos[jj].lon[jc][jr],kazm,filt_ar1[jc][jr]/filt_cnt[jc][jr],filt_dev[jc][jr],
		     grid_lat(grid_cell),grid_lon(grid_cell));
	      kazm_array=realloc(kazm_array,(ndata+1)*sizeof(struct mod_array));
	      kazm_array[ndata].coef=calloc(2*ngrid,sizeof(double));
	      los_data=realloc(los_data,(ndata+1)*sizeof(double));
	      los_err=realloc(los_err,(ndata+1)*sizeof(double));
	      kazm_array[ndata].coef[grid_cell]=sind(kazm);
	      kazm_array[ndata].coef[grid_cell+ngrid]=cosd(kazm);
	      los_data[ndata]=(-1.)*(filt_ar[jc][jr]+filt_ar1[jc][jr])/2;
	      los_err[ndata]=(filt_err[jc][jr]/filt_cnt[jc][jr])*(filt_err[jc][jr]/filt_cnt[jc][jr]);
	      ndata++;
	    }
	  }
	for( jr=0; jr<nbms; jr++)if(filt_ar[jr] != NULL)free(filt_ar[jr]); 
	free(filt_ar);
	for( jr=0; jr<nbms; jr++)if(filt_ar1[jr] != NULL)free(filt_ar1[jr]); 
	free(filt_ar1);
	for( jr=0; jr<nbms; jr++)if(filt_ar1[jr] != NULL)free(filt_dev[jr]); 
	free(filt_dev);
	for( jr=0; jr<nbms; jr++)if(filt_err[jr] != NULL)free(filt_err[jr]); 
	free(filt_err);
	for( jr=0; jr<nbms; jr++)if(filt_cnt[jr] != NULL)free(filt_cnt[jr]); 
	free(filt_cnt);
	for( jr=0; jr<nbms; jr++)if(good_ar[jr] != NULL)free(good_ar[jr]); 
	free(good_ar);
      }
    

      /* solve system */
      neqn=ngrid+nedge+ndata;
      if(smooth==1)neqn+=2*(ngrid-nedge);
      neqn+=(ngrid-nedge);

      fprintf(stderr,"NEQN=%d\n",neqn);
      coef=gsl_matrix_calloc((size_t)neqn,(size_t)ngrid*2);
      coefTrans=gsl_matrix_calloc((size_t)ngrid*2,(size_t)neqn);
      CD=gsl_matrix_calloc((size_t)neqn,(size_t)neqn);
      inv=gsl_matrix_calloc((size_t)neqn,(size_t)neqn);
      data=gsl_vector_calloc((size_t)neqn);


      fprintf(stderr,"***ALLOCATED***");
      jeqn=0;
      

      shf_fit_vel(grdfp,time,shf_vel,&shf_error);
      /*      ModelSigma=sqrt(2)/(shf_error*shf_error); */

      shfvf=fopen("shfvel_out","w");
      
      fprintf(shfvf,"%d\n",ngrid);
      for( jc=0; jc<ngrid; jc++)
	fprintf(shfvf,"%f %f %f %f %f\n",grid_lat(jc),grid_lon(jc),shf_vel[jc].vx,shf_vel[jc].vy,atan(shf_vel[jc].vx/shf_vel[jc].vy)/DTOR);
      
      fclose(shfvf);

	
      ModelSigma=model_scale*shf_error;
      for( jr=0; jr<nedge; jr++){
	for( jc=0; jc<2*ngrid; jc++)gsl_matrix_set(coef,jeqn,jc,edge_array[2*jr].coef[jc]);
	gsl_vector_set(data,jeqn,shf_vel[edge[jr]].vx);
	gsl_matrix_set(CD,jeqn,jeqn,ModelSigma);
	jeqn++;
	for( jc=0; jc<2*ngrid; jc++)gsl_matrix_set(coef,jeqn,jc,edge_array[2*jr+1].coef[jc]);
	gsl_vector_set(data,jeqn,shf_vel[edge[jr]].vy);
	gsl_matrix_set(CD,jeqn,jeqn,ModelSigma);
	jeqn++;
      }

      fprintf(stderr,"EDGES SET***");

      divscale=1.e5;
      for( jr=0; jr<ngrid-nedge; jr++){
	for( jc=0; jc<2*ngrid; jc++)gsl_matrix_set(coef,jeqn,jc,div_array[jr].coef[jc]*divscale);
	gsl_vector_set(data,jeqn,0);
	gsl_matrix_set(CD,jeqn,jeqn,10000.*MIN_ERR);
	jeqn++;
      }
      fprintf(stderr,"DIVERGENCE SET***");

      divscale=1.e5;
      for( jr=0; jr<ngrid-nedge; jr++){
	for( jc=0; jc<2*ngrid; jc++)gsl_matrix_set(coef,jeqn,jc,div_array1[jr].coef[jc]*divscale);
	gsl_vector_set(data,jeqn,0);
	gsl_matrix_set(CD,jeqn,jeqn,10000.*MIN_ERR);
	jeqn++;
      }
      fprintf(stderr,"DIVERGENCE 1 SET***");

      if(smooth==1){
	for( jr=0; jr<2*(ngrid-nedge); jr++){
	  for( jc=0; jc<2*ngrid; jc++)gsl_matrix_set(coef,jeqn,jc,smooth_array[jr].coef[jc]);
	  gsl_vector_set(data,jeqn,0);
	  gsl_matrix_set(CD,jeqn,jeqn,MIN_ERR);
	  jeqn++;
	}
	fprintf(stderr,"SMOOTHNESS SET***");
      }	

      for( jr=0; jr<ndata; jr++){
	for( jc=0; jc<2*ngrid; jc++)gsl_matrix_set(coef,jeqn,jc,kazm_array[jr].coef[jc]);
	gsl_vector_set(data,jeqn,los_data[jr]);
	gsl_matrix_set(CD,jeqn,jeqn,(MAX(los_err[jr],MIN_ERR)));
	jeqn++;
      }
      fprintf(stderr,"DATA SET***\n");

      /*      fp=fopen("coef_array","w");*/
      fprintf(fp,"%d %d %d\n",ngrid,nedge,neqn);
      for( jr=0; jr<neqn; jr++){
	colcount=0;
	for( jc=0; jc<2*ngrid; jc++)if( gsl_matrix_get(coef,jr,jc) !=0 )colcount++;
	fprintf(fp,"row %d %f %f %d\n",jr,gsl_vector_get(data,jr),gsl_matrix_get(CD,jr,jr),colcount);
	for( jc=0; jc<2*ngrid; jc++){
	  if( gsl_matrix_get(coef,jr,jc) !=0 ) fprintf(fp,"%d %f\n ",jc,gsl_matrix_get(coef,jr,jc));
	}
      }
      /*      fclose(fp); */

      GM=gsl_vector_calloc(neqn);
      for( jr=0; jr<neqn; jr++){
	temp=0.;
	for( jc=0; jc<ngrid; jc++ )temp+=gsl_matrix_get(coef,jr,jc)*shf_vel[jc].vx;
	for( jc=ngrid; jc<2*ngrid; jc++ )temp+=gsl_matrix_get(coef,jr,jc)*shf_vel[jc-ngrid].vy;
	gsl_vector_set(GM,jr,temp);
      }
      gsl_vector_sub(data,GM);
      gsl_vector_free(GM);

      t1=clock();
      gsl_matrix_transpose_memcpy(coefTrans,coef);
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,ModelSigma,coef,coefTrans,1.,CD);
      t2=clock();

      fprintf(stderr,"time for %d by %d gsl_blas_dgemm: %f\n",neqn,ngrid*2,((float)t2-(float)t1)/CLOCKS_PER_SEC);

      
      aa=malloc(neqn*sizeof(double*));
      for( jr=0; jr<neqn; jr++)aa[jr]=(double*)malloc(neqn*sizeof(double));
      
      for( jr=0; jr<neqn; jr++)for( jc=0; jc<neqn; jc++)aa[jr][jc]=gsl_matrix_get(CD,jr,jc);
      
      pinv_ar=malloc(neqn*sizeof(double*));
      for( jr=0; jr<neqn; jr++)pinv_ar[jr]=(double*)calloc(neqn,sizeof(double));

      if((s=pinv(neqn,neqn,aa,pinv_ar))!=-1){
	
	fprintf(stderr,"succesful return from pinv\n");
	for( jr=0; jr<neqn; jr++)for( jc=0; jc<neqn; jc++)gsl_matrix_set(inv,jr,jc,pinv_ar[jr][jc]);

	result=gsl_matrix_alloc(2*ngrid,neqn);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,ModelSigma,coefTrans,inv,0.,result);
	solution=gsl_vector_alloc(2*ngrid);
	gsl_blas_dgemv(CblasNoTrans,1.,result,data,0.,solution);
	gsl_matrix_free(result);
	
	
	printf("%d %d %d %f %f %f %f %d %f\n",0,0,0,0.,0.,0.,0.,0,0.);
	fprintf(vf,"%d %d %d %d %d %d\n",yr,mo,dy,hr,mt,(int)sc);
	fprintf(vf,"%f %f %f %f\n",0.,0.,0.,0.);
	fprintf(vf,"%d\n",ngrid);
	for( jc=0; jc<ngrid; jc++)
	  fprintf(vf,"%f %f %f %f %f %f %d\n",grid_lat(jc),grid_lon(jc),shf_vel[jc].vx+gsl_vector_get(solution,jc),shf_vel[jc].vy+gsl_vector_get(solution,jc+ngrid),shf_vel[jc].vx,shf_vel[jc].vy,grid_data_count[jc]);
	fflush(vf);	
	gsl_vector_free(solution);      

      }

      for( jr=0; jr<neqn; jr++)free(aa[jr]); 
      free(aa);
      for( jr=0; jr<neqn; jr++)free(pinv_ar[jr]); 
      free(pinv_ar);

      gsl_matrix_free(coef);
      gsl_matrix_free(coefTrans);
      gsl_matrix_free(CD);
      gsl_matrix_free(inv);
      gsl_vector_free(data);
      /*      kazm_array=realloc(kazm_array,sizeof(struct mod_array));
      //      los_data=realloc(los_data,sizeof(double));
      //      los_err=realloc(los_err,sizeof(double));
      */

      for( jr=0; jr<ndata; jr++)free(kazm_array[jr].coef);
      time+=avg_ival;
    }
  fclose(vf);
  fclose(fp); 
}
