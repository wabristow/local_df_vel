#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include "local_df_vel.h"


void cnv_to_str( int num, char *str)
{
  if(num < 10){sprintf(str,"0%d",num);} else{sprintf(str,"%d",num);}			
}

void find_fit_files(struct tm *t_s, struct tm *t_e, char* radar_list[], char* file_list[], int *nfiles)
{
  int yr,mo,dy;
  int i,nf=0;
  char yr_str[5], mo_str[3], dy_str[3];
  char *raid_path;
  char dir_path[256];
  DIR *dp;
  struct dirent *ep;

  raid_path=getenv("RAID_PATH");
  printf("%s\n",raid_path);

  int yr_st=t_s->tm_year+1900;
  int yr_nd=t_e->tm_year+1900;

  int mo_st=t_s->tm_mon+1;
  int mo_nd=t_e->tm_mon+1;

  int dy_st=t_s->tm_mday;
  int dy_nd=t_e->tm_mday;

  int yr_dy=yr_st*10000+mo_st*100+dy_st;
  int yr_dy_nd=yr_nd*10000+mo_nd*100+dy_nd;
  while( yr_dy <= yr_dy_nd ){
    yr=yr_dy/10000;
    mo=(yr_dy-10000*yr)/100;
    dy=(yr_dy-10000*yr-mo*100);

    CNV_TO_STR(yr,yr_str);
    CNV_TO_STR(mo,mo_str);
    CNV_TO_STR(dy,dy_str);

    sprintf(dir_path,"%s%s/%s.%s",raid_path,yr_str,mo_str,dy_str);
    printf("Directory is: %s\n",dir_path);
    dp=opendir(dir_path);
    while( (ep=readdir(dp)) )
      {
	i=-1;
	while( radar_list[++i] != NULL )
	  {
	    if( strstr(ep->d_name,".gz") != NULL) continue;
	    if( strstr(ep->d_name,".bz2") != NULL) continue;
	    if( strstr(ep->d_name,radar_list[i]) != NULL)
	      {
		nf++;
		file_list[nf]=malloc(sizeof(dir_path)+sizeof(ep->d_name));
		sprintf(file_list[nf],"%s/%s",dir_path,ep->d_name);
	      }
	  }
      }
    yr_dy++;
  }
  *nfiles=nf;
}
