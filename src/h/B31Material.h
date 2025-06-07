#pragma once

#include "Outputter.h"
#include "Material.h"

//!	Material class for B31 beam element
class CB31Material : public CMaterial
{
public:

    double A;   //!< Cross-sectional area
    double nu;   //！< Poisson ratio
    double G;   //!< Shear modulus
    double Iy;  //!< Moment of inertia about local y-axis
    double Iz;  //!< Moment of inertia about local z-axis
    double J;   //!< Torsional constant

public:

    //!	Read material data from stream Input
    virtual bool Read(ifstream& Input);

    //!	Write material data to Stream
    virtual void Write(COutputter& output);
};