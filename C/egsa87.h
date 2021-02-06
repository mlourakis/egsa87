#ifndef _EGSA87_H
#define _EGSA87_H

#ifdef __cplusplus
extern "C" {
#endif

extern void wgs84egsa87(double philam[2], double xy[2]);
extern void egsa87wgs84(double xy[2], double philam[2]);

#ifdef __cplusplus
}
#endif

#endif /* _EGSA87_H */

