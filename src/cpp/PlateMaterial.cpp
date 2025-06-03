#include "PlateMaterial.h"
#include <iomanip>

using namespace std;

bool CPlateMaterial::Read(ifstream& Input)
{
    Input >> nset >> E >> nu >> thickness;
    return true;
}

void CPlateMaterial::Write(COutputter& output)
{
    output << setw(16) << E
           << setw(16) << nu
           << setw(16) << thickness << endl;
}

Eigen::Matrix3d CPlateMaterial::GetFlexuralMatrix() const {
    double coeff = E * thickness * thickness * thickness / (12.0 * (1.0 - nu * nu));
    Eigen::Matrix3d D;
    D << coeff,     coeff*nu,  0.0,
         coeff*nu,  coeff,     0.0,
         0.0,       0.0,       coeff*(1.0-nu)/2.0;
    return D;
}


Eigen::Matrix3d CPlateMaterial::GetstressMatrix() const {
    double coeff = E / (1.0 - nu * nu);
    Eigen::Matrix3d D;
    D << coeff,     coeff*nu,  0.0,
         coeff*nu,  coeff,     0.0,
         0.0,       0.0,       coeff*(1.0-nu)/2.0;
    return D;
}