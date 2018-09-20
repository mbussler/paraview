#include "densityblob.h"

//---Blob---------------------------------------------------------------------
double DensityBlob::GetDensity(double x, double y, double z) {
    double distsqr = (m_cx-x)*(m_cx-x)+(m_cy-y)*(m_cy-y)+(m_cz-z)*(m_cz-z);
    if (sqrt(distsqr) < m_r )
	return ( pow(m_s*(1-(distsqr/m_r/m_r)) , 2) );
    else
	return 0.0;
    //return sqrt(distsqr);
}
