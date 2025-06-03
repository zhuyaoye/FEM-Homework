#pragma once

#include "Outputter.h"
#include "Material.h"

//! Material class for C3D8R 8-node linear brick element with reduced integration
class CC3D8RMaterial : public CMaterial
{
public:
    double E;    // Young's modulus
    double nu;   // Poisson's ratio

public:
    //! Read material data from stream Input
    virtual bool Read(ifstream& Input);

    //! Write material data to Stream
    virtual void Write(COutputter& output);
};