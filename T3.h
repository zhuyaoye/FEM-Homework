#pragma once

#include "Element.h"


class CT3 : public CElement
{
public:

	//!	Constructor
	CT3();

	//!	Desconstructor
	~CT3();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output);

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

private:
   
    void ShapeFunction(double* N, double* dNdx, double* dNdy); // shape function
	double CalculateArea();
	
};