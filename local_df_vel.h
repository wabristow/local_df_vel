#include "cnvgrid.h"

#define PATH_LEN 512
#define F_NAME_LEN 64
#define MXFILES 1000
#define NRAD 33
#define MAX_RANGE 1000
#define MAX_BEAM 25

#define CNV_TO_STR(num,str) if(num < 10){sprintf(str,"0%d",num);} else{sprintf(str,"%d",num);}
#define IS_LEAPYEAR(year) ((year%4)&&(!(year%100))||(year%400))

typedef struct cell{
  double lat[4];
  double lon[4];
  int edge;
}CELL;

typedef struct boundary{
  double lat;
  double lon;
}B_POINT;

typedef struct Neighbor{
  int c[4];
}NEIGHBOR;

typedef struct mod_array{
  double* coef;
}M_ARRAY;

typedef struct time{
  int yr;
  int mo;
  int dy;
  int hr;
  int mt;
  int sc;
  long us;
} TIME;

typedef struct finfo{
  char dir_path[PATH_LEN];
  char fname[F_NAME_LEN];
}FILE_INFO;

typedef struct RadarPos{
  double lat[MAX_BEAM][MAX_RANGE];
  double lon[MAX_BEAM][MAX_RANGE];
  double kazm[MAX_BEAM][MAX_RANGE];
}RPOS;  

typedef struct RadarMap{
  int cell[MAX_BEAM][MAX_RANGE];
}RMAP;  

typedef struct PotRec{
  float pot;
  int count;
  float glat;
  float glon;
  float mlat;
  float mlon;
  float mlt;
}POTEN;

typedef struct ShfVel{
  double vx;
  double vy;
}SHFV;

void parse_instructions();
/* void find_fit_files(struct tm *,struct tm *,char *[], char *[], int *); */
struct tm *parse_date_str(long);

int make_local_pgrid(struct cell * , int, struct CnvGrid * );

typedef struct VDSTRS{
  int yr;
  int mo;
  int dy;
  int hr;
  int mt;
  float sc;
  double Bx,By,Bz;
  double delay;
  int npts;
  double* lat;
  double* lon;
  double* vx;
  double* vy;
  double* ex;
  double* ey;
  double* vaz;
  double* vmag;
  double* mvaz;
  double* mvmag;
  int* dcount;
}VDSTR;

