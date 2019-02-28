
#include "astalg.h"
#include "rtime.h"

struct TPOS {
  double tlat;
  double tlon;
};

int sub_sphazm(double Alon,double Alat,double Clon,double Clat,double *azm,double *range);
int sub_sphcal(double Alon, double Alat, double azm, double range, double *Clon, double *Clat);

