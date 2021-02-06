#include <stdio.h>
#include <math.h>

#include "egsa87.h"


int main()
{
// https://en.wikipedia.org/wiki/Morosini_Fountain
double phlam[2]={35.3391370, 25.1331727};
double xy[2], phlam2[2];

  phlam[0]*=M_PI/180.0;
  phlam[1]*=M_PI/180.0;
  wgs84egsa87(phlam, xy);
  printf("Point in EGSA87: %lf %lf\n", xy[0], xy[1]);

  // back to WGS84
  egsa87wgs84(xy, phlam2);
  phlam2[0]*=180.0/M_PI;
  phlam2[1]*=180.0/M_PI;
  printf("Point in WGS84: %lf %lf\n", phlam2[0], phlam2[1]);
}
