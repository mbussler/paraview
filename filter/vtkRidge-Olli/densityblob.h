#ifndef _DENSITYBLOB_H_
#define _DENSITYBLOB_H_

#include <math.h>

class DensityBlob {
 public:
    inline void SetCenter(double x, double y, double z) {m_cx = x; m_cy = y; m_cz = z;}
    inline void SetRadius(double radius) {m_r = radius;}
    inline void SetStrenght(double strenght) {m_s = strenght;}
    double GetDensity(double,double,double);
 private:
    double m_cx,m_cy,m_cz;
    double m_r,m_s;
};

#endif //_DENSITYBLOB_H_
