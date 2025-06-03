#pragma once

#include "Element.h"

//! C3D8R 8-node linear brick element with reduced integration
class CC3D8R : public CElement
{
public:
    //! Constructor
    CC3D8R();

    //! Destructor
    ~CC3D8R();

    //! Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

    //! Write element data to stream
    virtual void Write(COutputter& output);

    //! Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

    //! Calculate element stress
    virtual void ElementStress(double* stress, double* Displacement);

private:
    //! Material properties
    double E;       // Young's modulus
    double nu;      // Poisson's ratio
    
    //! Shape function and derivatives
    void ShapeFunction(double xi, double eta, double zeta, double N[8]);
    void ShapeFunctionDerivatives(double xi, double eta, double zeta, double dN[8][3]);
    
    //! Calculate Jacobian matrix and its inverse
    void Jacobian(double xi, double eta, double zeta, double J[3][3], double& detJ);
    
    //! Calculate B matrix (strain-displacement matrix)
    void BMatrix(double xi, double eta, double zeta, double B[6][24]);
    
    //! Calculate D matrix (constitutive matrix)
    void DMatrix(double D[6][6]);
    
    //! Gauss integration points and weights (1-point reduced integration)
    static const int NumGaussPoints = 1;
    static const double GaussPoints[NumGaussPoints];
    static const double GaussWeights[NumGaussPoints];
};