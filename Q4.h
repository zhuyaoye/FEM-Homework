#pragma once

#include "Element.h"


class CQ4: public CElement
{
public:

	//!	Constructor
	CQ4();

	//!	Desconstructor
	~CQ4();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output);

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

private:

    void ShapeFunction(double xi, double eta, double* N, double* dNdxi, double* dNdeta); // �κ���
    void Jacobian(double xi, double eta, double& detJ, double* dNdx, double* dNdy);     // �ſɱȾ���
    void GaussIntegration(double* Matrix);                                              // ��˹����
    void ConstructBMatrix(double B[3][8], const double* dNdx, const double* dNdy);      // ����B����
};