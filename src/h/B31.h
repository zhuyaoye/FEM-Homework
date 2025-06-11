#pragma once

#include "Element.h"
#include "B31Material.h"

//! B31 Beam element class
class CB31 : public CElement
{
public:

//!	Constructor
	CB31();

//!	Desconstructor
	~CB31();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

private:

//!	Calculate transformation matrix
	void CalculateTransformationMatrix(double T[12][12]);

//! Length of the beam
	double Length_;

//! Direction vector of the beam
	double Direction_[3];

//! Referance vector of the beam (using to caculate the local y, ordinary global z)
	double UpVector_[3];
//! Material properties
	double E, nu, G, A, Iy, Iz, J;
};