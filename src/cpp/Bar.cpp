/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Bar.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CBar::CBar()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CBar::~CBar()
{
}

//	Read element data from stream Input
bool CBar::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
    ElementMaterial_ = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CBar::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBar::ElementStiffness(double* Matrix)
{
    clear(Matrix, 144);  // 12x12 full DOF matrix

    // Compute direction vector
    double DX[3];
    for (unsigned int i = 0; i < 3; i++)
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

    double DX2[6];
    DX2[0] = DX[0] * DX[0];
    DX2[1] = DX[1] * DX[1];
    DX2[2] = DX[2] * DX[2];
    DX2[3] = DX[0] * DX[1];
    DX2[4] = DX[1] * DX[2];
    DX2[5] = DX[0] * DX[2];

    double L2 = DX2[0] + DX2[1] + DX2[2];
    double L = sqrt(L2);

    CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);
    double k = material_->E * material_->Area / (L * L2);

    // Original 6x6 stiffness matrix
    double K_small[6][6] = {0};

    K_small[0][0] =  k * DX2[0];
    K_small[0][1] =  k * DX2[3];
    K_small[0][2] =  k * DX2[5];
    K_small[0][3] = -k * DX2[0];
    K_small[0][4] = -k * DX2[3];
    K_small[0][5] = -k * DX2[5];

    K_small[1][1] =  k * DX2[1];
    K_small[1][2] =  k * DX2[4];
    K_small[1][3] = -k * DX2[3];
    K_small[1][4] = -k * DX2[1];
    K_small[1][5] = -k * DX2[4];

    K_small[2][2] =  k * DX2[2];
    K_small[2][3] = -k * DX2[5];
    K_small[2][4] = -k * DX2[4];
    K_small[2][5] = -k * DX2[2];

    K_small[3][3] =  k * DX2[0];
    K_small[3][4] =  k * DX2[3];
    K_small[3][5] =  k * DX2[5];

    K_small[4][4] =  k * DX2[1];
    K_small[4][5] =  k * DX2[4];

    K_small[5][5] =  k * DX2[2];

    // Mirror lower triangle to upper triangle
    for (int i = 0; i < 6; ++i)
        for (int j = i + 1; j < 6; ++j)
            K_small[j][i] = K_small[i][j];

    // Full 12x12 matrix (initialized to zero)
    double K_full[12][12] = {0};

    // Map 6x6 small K into full 12x12 K_full
    // Node 1: DOF 0~2 (UX, UY, UZ), Node 2: DOF 6~8 (UX, UY, UZ)
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            K_full[i][j]       = K_small[i][j];       // K11
            K_full[i][j + 6]   = K_small[i][j + 3];   // K12
            K_full[i + 6][j]   = K_small[i + 3][j];   // K21
            K_full[i + 6][j + 6] = K_small[i + 3][j + 3]; // K22
        }
    }
	int index = 0;
    const int n = 12;
    for (int k = 0; k < n; ++k)
    {
        Matrix[index++] = K_full[k][k];
        for (int i = k - 1; i >= 0; --i)
        {
            Matrix[index++] = K_full[i][k];
        }
    }
}

//	Calculate element stress 
void CBar::ElementStress(double* stress, double* Displacement)
{
	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i+3] = -S[i];
	}
	
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += S[i] * Displacement[LocationMatrix_[i]-1];
	}
}
