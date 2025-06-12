#include "C3D8R.h"
#include "C3D8RMaterial.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// Gaussian integration points (reduced integration, 1 point)
const double CC3D8R::GaussPoints[NumGaussPoints] = {0.0};
const double CC3D8R::GaussWeights[NumGaussPoints] = {8.0};

// Constructor
CC3D8R::CC3D8R()
{
    NEN_ = 8;    // 8 nodes
    ND_ = 24;    // 3 DOF per node × 8
    nodes_ = new CNode*[NEN_];
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}

// Destructor
CC3D8R::~CC3D8R()
{
    delete[] nodes_;
    delete[] LocationMatrix_;
}

// Read element data from input stream
bool CC3D8R::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;
    Input >> MSet;  // Material set number
    
    // Read 8 node numbers
    for(int i=0; i<8; i++) {
        unsigned int N;
        Input >> N;
        nodes_[i] = &NodeList[N-1];
    }
    
    ElementMaterial_ = dynamic_cast<CC3D8RMaterial*>(MaterialSets + MSet - 1);
    return true;
}

// Write element information to output stream
void CC3D8R::Write(COutputter& output)
{
    output << setw(9) << ElementMaterial_->nset;
    for(int i=0; i<8; i++) {
        output << setw(9) << nodes_[i]->NodeNumber;
    }
    output << endl;
}

// Compute element stiffness matrix (24×24)
void CC3D8R::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());
    
    CC3D8RMaterial* material = dynamic_cast<CC3D8RMaterial*>(ElementMaterial_);
    double D[6][6];
    DMatrix(D);  // Compute constitutive matrix
    
    // Reduced integration - single Gauss point
    for(int gp=0; gp<NumGaussPoints; gp++) {
        double xi = GaussPoints[gp];
        double eta = GaussPoints[gp];
        double zeta = GaussPoints[gp];
        double weight = GaussWeights[gp];
        
        double B[6][24];
        BMatrix(xi, eta, zeta, B);  // Compute B matrix
        
        double J[3][3], detJ;
        Jacobian(xi, eta, zeta, J, detJ);  // Compute Jacobian matrix
        
        // K += B^T * D * B * detJ * weight
        for(int i=0; i<24; i++) {
            for(int j=0; j<24; j++) {
                double sum = 0.0;
                for(int k=0; k<6; k++) {
                    for(int l=0; l<6; l++) {
                        sum += B[k][i] * D[k][l] * B[l][j];
                    }
                }
                Matrix[i*24 + j] += sum * detJ * weight;
            }
        }
    }
}

// Calculate element stress
void CC3D8R::ElementStress(double* stress, double* Displacement)
{
    // Use center point for stress calculation (reduced integration)
    double xi = 0.0, eta = 0.0, zeta = 0.0;
    
    double B[6][24];
    BMatrix(xi, eta, zeta, B);
    
    // Strain = B * U
    double strain[6] = {0.0};
    for(int i=0; i<6; i++) {
        for(int j=0; j<24; j++) {
            strain[i] += B[i][j] * Displacement[LocationMatrix_[j]-1];
        }
    }
    
    // Stress = D * strain
    double D[6][6];
    DMatrix(D);
    
    for(int i=0; i<6; i++) {
        stress[i] = 0.0;
        for(int j=0; j<6; j++) {
            stress[i] += D[i][j] * strain[j];
        }
    }
}

// Compute shape functions at given natural coordinates
void CC3D8R::ShapeFunction(double xi, double eta, double zeta, double N[8])
{
    N[0] = 0.125*(1-xi)*(1-eta)*(1-zeta);
    N[1] = 0.125*(1+xi)*(1-eta)*(1-zeta);
    N[2] = 0.125*(1+xi)*(1+eta)*(1-zeta);
    N[3] = 0.125*(1-xi)*(1+eta)*(1-zeta);
    N[4] = 0.125*(1-xi)*(1-eta)*(1+zeta);
    N[5] = 0.125*(1+xi)*(1-eta)*(1+zeta);
    N[6] = 0.125*(1+xi)*(1+eta)*(1+zeta);
    N[7] = 0.125*(1-xi)*(1+eta)*(1+zeta);
}

// Compute derivatives of shape functions with respect to natural coordinates
void CC3D8R::ShapeFunctionDerivatives(double xi, double eta, double zeta, double dN[8][3])
{
    // dN/dxi
    dN[0][0] = -0.125*(1-eta)*(1-zeta);
    dN[1][0] =  0.125*(1-eta)*(1-zeta);
    dN[2][0] =  0.125*(1+eta)*(1-zeta);
    dN[3][0] = -0.125*(1+eta)*(1-zeta);
    dN[4][0] = -0.125*(1-eta)*(1+zeta);
    dN[5][0] =  0.125*(1-eta)*(1+zeta);
    dN[6][0] =  0.125*(1+eta)*(1+zeta);
    dN[7][0] = -0.125*(1+eta)*(1+zeta);
    
    // dN/deta
    dN[0][1] = -0.125*(1-xi)*(1-zeta);
    dN[1][1] = -0.125*(1+xi)*(1-zeta);
    dN[2][1] =  0.125*(1+xi)*(1-zeta);
    dN[3][1] =  0.125*(1-xi)*(1-zeta);
    dN[4][1] = -0.125*(1-xi)*(1+zeta);
    dN[5][1] = -0.125*(1+xi)*(1+zeta);
    dN[6][1] =  0.125*(1+xi)*(1+zeta);
    dN[7][1] =  0.125*(1-xi)*(1+zeta);
    
    // dN/dzeta
    dN[0][2] = -0.125*(1-xi)*(1-eta);
    dN[1][2] = -0.125*(1+xi)*(1-eta);
    dN[2][2] = -0.125*(1+xi)*(1+eta);
    dN[3][2] = -0.125*(1-xi)*(1+eta);
    dN[4][2] =  0.125*(1-xi)*(1-eta);
    dN[5][2] =  0.125*(1+xi)*(1-eta);
    dN[6][2] =  0.125*(1+xi)*(1+eta);
    dN[7][2] =  0.125*(1-xi)*(1+eta);
}

// Compute Jacobian matrix and its determinant
void CC3D8R::Jacobian(double xi, double eta, double zeta, double J[3][3], double& detJ)
{
    double dN[8][3];
    ShapeFunctionDerivatives(xi, eta, zeta, dN);
    
    // Initialize Jacobian matrix
    for(int i=0; i<3; i++) {
        for(int j=0; j<3; j++) {
            J[i][j] = 0.0;
        }
    }
    
    // Compute Jacobian matrix components
    for(int i=0; i<8; i++) {
        double x = nodes_[i]->XYZ[0];
        double y = nodes_[i]->XYZ[1];
        double z = nodes_[i]->XYZ[2];
        
        for(int j=0; j<3; j++) {
            J[0][j] += x * dN[i][j];
            J[1][j] += y * dN[i][j];
            J[2][j] += z * dN[i][j];
        }
    }
    
    // Compute determinant of Jacobian
    detJ = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) -
           J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0]) +
           J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
}

// Compute B matrix (strain-displacement matrix)
void CC3D8R::BMatrix(double xi, double eta, double zeta, double B[6][24])
{
    double dN[8][3];
    ShapeFunctionDerivatives(xi, eta, zeta, dN);
    
    double J[3][3], detJ;
    Jacobian(xi, eta, zeta, J, detJ);
    
    // Compute inverse of Jacobian matrix
    double invJ[3][3];
    double invDetJ = 1.0 / detJ;
    
    invJ[0][0] =  (J[1][1]*J[2][2] - J[1][2]*J[2][1]) * invDetJ;
    invJ[0][1] = -(J[0][1]*J[2][2] - J[0][2]*J[2][1]) * invDetJ;
    invJ[0][2] =  (J[0][1]*J[1][2] - J[0][2]*J[1][1]) * invDetJ;
    invJ[1][0] = -(J[1][0]*J[2][2] - J[1][2]*J[2][0]) * invDetJ;
    invJ[1][1] =  (J[0][0]*J[2][2] - J[0][2]*J[2][0]) * invDetJ;
    invJ[1][2] = -(J[0][0]*J[1][2] - J[0][2]*J[1][0]) * invDetJ;
    invJ[2][0] =  (J[1][0]*J[2][1] - J[1][1]*J[2][0]) * invDetJ;
    invJ[2][1] = -(J[0][0]*J[2][1] - J[0][1]*J[2][0]) * invDetJ;
    invJ[2][2] =  (J[0][0]*J[1][1] - J[0][1]*J[1][0]) * invDetJ;
    
    // Initialize B matrix
    for(int i=0; i<6; i++) {
        for(int j=0; j<24; j++) {
            B[i][j] = 0.0;
        }
    }
    
    // Compute B matrix components
    for(int i=0; i<8; i++) {
        // Derivatives of shape functions with respect to global coordinates
        double dNdx = dN[i][0]*invJ[0][0] + dN[i][1]*invJ[0][1] + dN[i][2]*invJ[0][2];
        double dNdy = dN[i][0]*invJ[1][0] + dN[i][1]*invJ[1][1] + dN[i][2]*invJ[1][2];
        double dNdz = dN[i][0]*invJ[2][0] + dN[i][1]*invJ[2][1] + dN[i][2]*invJ[2][2];
        
        int col = i*3;
        
        // Strain-displacement relations
        B[0][col]   = dNdx;  // εxx
        B[1][col+1] = dNdy;  // εyy
        B[2][col+2] = dNdz;  // εzz
        
        B[3][col]   = dNdy;  // γxy
        B[3][col+1] = dNdx;
        
        B[4][col+1] = dNdz;  // γyz
        B[4][col+2] = dNdy;
        
        B[5][col]   = dNdz;  // γzx
        B[5][col+2] = dNdx;
    }
}

// Compute D matrix (constitutive matrix for isotropic linear elasticity)
void CC3D8R::DMatrix(double D[6][6])
{
    CC3D8RMaterial* material = dynamic_cast<CC3D8RMaterial*>(ElementMaterial_);
    double E = material->E;
    double nu = material->nu;
    
    double lambda = E*nu/((1+nu)*(1-2*nu));
    double mu = E/(2*(1+nu));
    
    // Initialize D matrix
    for(int i=0; i<6; i++) {
        for(int j=0; j<6; j++) {
            D[i][j] = 0.0;
        }
    }
    
    // Fill D matrix components (isotropic linear elasticity)
    D[0][0] = D[1][1] = D[2][2] = lambda + 2*mu;
    D[0][1] = D[1][0] = D[0][2] = D[2][0] = D[1][2] = D[2][1] = lambda;
    D[3][3] = D[4][4] = D[5][5] = mu;
}
