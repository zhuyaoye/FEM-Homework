#pragma once

#include "Material.h"
#include "Outputter.h"
#include <Eigen/Dense>
class CPlateMaterial : public CMaterial {
public:
    // Read material data from stream
    virtual bool Read(std::ifstream& Input);

    // Write material data to output stream
    virtual void Write(COutputter& output);

    // Flexural matrix D
    Eigen::Matrix3d GetFlexuralMatrix() const;

    // 平面应力D
    Eigen::Matrix3d GetstressMatrix() const;

    // Getters
    double GetE() const { return E; }
    double GetNu() const { return nu; }
    double GetThickness() const { return thickness; }

private:
    double E;         // Young's modulus
    double nu;        // Poisson's ratio
    double thickness; // Plate thickness
};