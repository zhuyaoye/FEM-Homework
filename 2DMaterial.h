
#pragma once

#include "Outputter.h"
#include "Material.h"



// ! Material class for T3 element
class C2DMaterial : public CMaterial
{
public:
   
    double nu;   // Poisson ratio
    double t;    // thick
    int type;    // 0-plane  stress 1-plane strain


public:

    //!	Read material data from stream Input
    virtual bool Read(ifstream& Input);

    //!	Write material data to Stream
    virtual void Write(COutputter& output);
};