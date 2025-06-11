#include "B31.h"
#include "B31Material.h"
#include <iostream>
#include <iomanip>
#include <cmath>

#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// Constructor
CB31::CB31()
{
    NEN_ = 2;    // 2 nodes
    ND_ = 12;    // 6 DOF per node × 2
    nodes_ = new CNode*[NEN_];
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}

// Destructor
CB31::~CB31()
{
}

bool CB31::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int N1, N2, MSet;
    std::string line;
    if (!std::getline(Input, line))
        return false;  // 文件结束或读取失败
    
    std::istringstream iss(line);
    
    if (!(iss >> N1 >> N2 >> MSet))
        return false;  // 读取节点和材料集失败
    
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];

    ElementMaterial_ = dynamic_cast<CB31Material*>(MaterialSets) + MSet - 1;

    // 检查是否有额外数据用于up向量
    double upx, upy, upz;
    if (iss >> upx >> upy >> upz) {
        UpVector_[0] = upx;
        UpVector_[1] = upy;
        UpVector_[2] = upz;
    } else {
        // 没有提供up向量，使用默认值
        UpVector_[0] = 0.0;
        UpVector_[1] = 0.0;
        UpVector_[2] = 1.0;
    }
    
    return true;
}

// Write element info
void CB31::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9)  << nodes_[1]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

// Compute element stiffness matrix (12×12)
void CB31::ElementStiffness(double* Matrix)
{
      clear(Matrix, SizeOfStiffnessMatrix());

    //	Calculate bar length
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double L2 = DX[0] *  DX[0] + DX[1] * DX[1] + DX[2] * DX[2];
	double L = sqrt(L2);
    double L3 = L * L2;

    //	Calculate element stiffness matrix

	CB31Material* material_ = dynamic_cast<CB31Material*>(ElementMaterial_);	// Pointer to material of the element

    // Calculate the coordinate transformation matrix
    double Z_direction[3];
    for(unsigned i =0; i < 3; i++)
        Z_direction[i] = UpVector_[i];
    
    double L_Z = sqrt(Z_direction[0] * Z_direction[0] + Z_direction[1] * Z_direction[1] + + Z_direction[2] * Z_direction[2]);

    for(unsigned i =0; i < 3; i++)
        Z_direction[i] = Z_direction[i]/L_Z;
    
    double L00 = DX[0]/L;
    double L01 = DX[1]/L;
    double L02 = DX[2]/L;
    double L20 = Z_direction[0];
    double L21 = Z_direction[1];
    double L22 = Z_direction[2];

    double DOT;
    DOT = L00*L20 + L01*L21 + L02*L22;

    if(DOT != 0)
    {
        cerr << "The direction of top surface or the coordinate of certain node is wrong" << endl;
        exit(EXIT_FAILURE);
    }

    double L10 = L21 * L02 - L22 * L01;
    double L11 = L22 * L00 - L20 * L02;
    double L12 = L20 * L01 - L21 * L00;

    double Area = material_->A; //Area of the beam's section
    double k1 = material_->E * Area / L; // Coefficient of axial tension and compression

    double G = material_->E / 2 / (1 + material_->nu); //  Calculate the Shear Modulus
    double J = material_->J;  // Polar moment of inertia
    double k2 = G * J / L;  // Coefficient of torsional stiffness

    double I_y = material_->Iy; // Calculate the moment of inertia of plate xOz
    double I_z = material_->Iz; // Calculate the moment of inertia of plate xOy
    double k3 = material_->E * I_y / L3; // Coefficient of the beam bending in the y-direction
    double k4 = material_->E * I_z / L3; // Coefficient of the beam bending in the z-direction


    // 列优先顺序，仅上半部分（含对角线）
    Matrix[0] = k1 * pow(L00, 2) + 12 * k4 * pow(L10, 2) + 12 * k3 * pow(L20, 2);  // 列1: 对角元

    Matrix[1] = k1 * pow(L01, 2) + 12 * k4 * pow(L11, 2) + 12 * k3 * pow(L21, 2);  // 列2: 对角元
    Matrix[2] = L00 * L01 * k1 + 12 * L10 * L11 * k4 + 12 * L20 * L21 * k3;  // 列2: 第1行

    Matrix[3] = k1 * pow(L02, 2) + 12 * k4 * pow(L12, 2) + 12 * k3 * pow(L22, 2);  // 列3: 对角元
    Matrix[4] = L01 * L02 * k1 + 12 * L11 * L12 * k4 + 12 * L21 * L22 * k3;  // 列3: 第2行
    Matrix[5] = L00 * L02 * k1 + 12 * L10 * L12 * k4 + 12 * L20 * L22 * k3;  // 列3: 第1行

    Matrix[6] = 4 * k3 * pow(L, 2) * pow(L10, 2) + 4 * k4 * pow(L, 2) * pow(L20, 2) + k2 * pow(L00, 2);  // 列4: 对角元
    Matrix[7] = 6 * L * L12 * L20 * k4 - 6 * L * L10 * L22 * k3;  // 列4: 第3行
    Matrix[8] = 6 * L * L11 * L20 * k4 - 6 * L * L10 * L21 * k3;  // 列4: 第2行
    Matrix[9] = 6 * L * L10 * L20 * k4 - 6 * L * L10 * L20 * k3;  // 列4: 第1行

    Matrix[10] = 4 * k3 * pow(L, 2) * pow(L11, 2) + 4 * k4 * pow(L, 2) * pow(L21, 2) + k2 * pow(L01, 2);  // 列5: 对角元
    Matrix[11] = L00 * L01 * k2 + 4 * pow(L, 2) * L10 * L11 * k3 + 4 * pow(L, 2) * L20 * L21 * k4;  // 列5: 第4行
    Matrix[12] = 6 * L * L12 * L21 * k4 - 6 * L * L11 * L22 * k3;  // 列5: 第3行
    Matrix[13] = 6 * L * L11 * L21 * k4 - 6 * L * L11 * L21 * k3;  // 列5: 第2行
    Matrix[14] = 6 * L * L10 * L21 * k4 - 6 * L * L11 * L20 * k3;  // 列5: 第1行

    Matrix[15] = 4 * k3 * pow(L, 2) * pow(L12, 2) + 4 * k4 * pow(L, 2) * pow(L22, 2) + k2 * pow(L02, 2);  // 列6: 对角元
    Matrix[16] = L01 * L02 * k2 + 4 * pow(L, 2) * L11 * L12 * k3 + 4 * pow(L, 2) * L21 * L22 * k4;  // 列6: 第5行
    Matrix[17] = L00 * L02 * k2 + 4 * pow(L, 2) * L10 * L12 * k3 + 4 * pow(L, 2) * L20 * L22 * k4;  // 列6: 第4行
    Matrix[18] = 6 * L * L12 * L22 * k4 - 6 * L * L12 * L22 * k3;  // 列6: 第3行
    Matrix[19] = 6 * L * L11 * L22 * k4 - 6 * L * L12 * L21 * k3;  // 列6: 第2行
    Matrix[20] = 6 * L * L10 * L22 * k4 - 6 * L * L12 * L20 * k3;  // 列6: 第1行

    Matrix[21] = k1 * pow(L00, 2) + 12 * k4 * pow(L10, 2) + 12 * k3 * pow(L20, 2);  // 列7: 对角元
    Matrix[22] = 6 * L * L12 * L20 * k3 - 6 * L * L10 * L22 * k4;  // 列7: 第6行
    Matrix[23] = 6 * L * L11 * L20 * k3 - 6 * L * L10 * L21 * k4;  // 列7: 第5行
    Matrix[24] = 6 * L * L10 * L20 * k3 - 6 * L * L10 * L20 * k4;  // 列7: 第4行
    Matrix[25] = -L00 * L02 * k1 - 12 * L10 * L12 * k4 - 12 * L20 * L22 * k3;  // 列7: 第3行
    Matrix[26] = -L00 * L01 * k1 - 12 * L10 * L11 * k4 - 12 * L20 * L21 * k3;  // 列7: 第2行
    Matrix[27] = -k1 * pow(L00, 2) - 12 * k4 * pow(L10, 2) - 12 * k3 * pow(L20, 2);  // 列7: 第1行

    Matrix[28] = k1 * pow(L01, 2) + 12 * k4 * pow(L11, 2) + 12 * k3 * pow(L21, 2);  // 列8: 对角元
    Matrix[29] = L00 * L01 * k1 + 12 * L10 * L11 * k4 + 12 * L20 * L21 * k3;  // 列8: 第7行
    Matrix[30] = 6 * L * L12 * L21 * k3 - 6 * L * L11 * L22 * k4;  // 列8: 第6行
    Matrix[31] = 6 * L * L11 * L21 * k3 - 6 * L * L11 * L21 * k4;  // 列8: 第5行
    Matrix[32] = 6 * L * L10 * L21 * k3 - 6 * L * L11 * L20 * k4;  // 列8: 第4行
    Matrix[33] = -L01 * L02 * k1 - 12 * L11 * L12 * k4 - 12 * L21 * L22 * k3;  // 列8: 第3行
    Matrix[34] = -k1 * pow(L01, 2) - 12 * k4 * pow(L11, 2) - 12 * k3 * pow(L21, 2);  // 列8: 第2行
    Matrix[35] = -L00 * L01 * k1 - 12 * L10 * L11 * k4 - 12 * L20 * L21 * k3;  // 列8: 第1行

    Matrix[36] = k1 * pow(L02, 2) + 12 * k4 * pow(L12, 2) + 12 * k3 * pow(L22, 2);  // 列9: 对角元
    Matrix[37] = L01 * L02 * k1 + 12 * L11 * L12 * k4 + 12 * L21 * L22 * k3;  // 列9: 第8行
    Matrix[38] = L00 * L02 * k1 + 12 * L10 * L12 * k4 + 12 * L20 * L22 * k3;  // 列9: 第7行
    Matrix[39] = 6 * L * L12 * L22 * k3 - 6 * L * L12 * L22 * k4;  // 列9: 第6行
    Matrix[40] = 6 * L * L11 * L22 * k3 - 6 * L * L12 * L21 * k4;  // 列9: 第5行
    Matrix[41] = 6 * L * L10 * L22 * k3 - 6 * L * L12 * L20 * k4;  // 列9: 第4行
    Matrix[42] = -k1 * pow(L02, 2) - 12 * k4 * pow(L12, 2) - 12 * k3 * pow(L22, 2);  // 列9: 第3行
    Matrix[43] = -L01 * L02 * k1 - 12 * L11 * L12 * k4 - 12 * L21 * L22 * k3;  // 列9: 第2行
    Matrix[44] = -L00 * L02 * k1 - 12 * L10 * L12 * k4 - 12 * L20 * L22 * k3;  // 列9: 第1行

    Matrix[45] = 4 * k3 * pow(L, 2) * pow(L10, 2) + 4 * k4 * pow(L, 2) * pow(L20, 2) + k2 * pow(L00, 2);  // 列10: 对角元
    Matrix[46] = 6 * L * L10 * L22 * k3 - 6 * L * L12 * L20 * k4;  // 列10: 第9行
    Matrix[47] = 6 * L * L10 * L21 * k3 - 6 * L * L11 * L20 * k4;  // 列10: 第8行
    Matrix[48] = 6 * L * L10 * L20 * k3 - 6 * L * L10 * L20 * k4;  // 列10: 第7行
    Matrix[49] = 2 * pow(L, 2) * L10 * L12 * k3 - L00 * L02 * k2 + 2 * pow(L, 2) * L20 * L22 * k4;  // 列10: 第6行
    Matrix[50] = 2 * pow(L, 2) * L10 * L11 * k3 - L00 * L01 * k2 + 2 * pow(L, 2) * L20 * L21 * k4;  // 列10: 第5行
    Matrix[51] = 2 * k3 * pow(L, 2) * pow(L10, 2) + 2 * k4 * pow(L, 2) * pow(L20, 2) - k2 * pow(L00, 2);  // 列10: 第4行
    Matrix[52] = 6 * L * L12 * L20 * k4 - 6 * L * L10 * L22 * k3;  // 列10: 第3行
    Matrix[53] = 6 * L * L11 * L20 * k4 - 6 * L * L10 * L21 * k3;  // 列10: 第2行
    Matrix[54] = 6 * L * L10 * L20 * k4 - 6 * L * L10 * L20 * k3;  // 列10: 第1行

    Matrix[55] = 4 * k3 * pow(L, 2) * pow(L11, 2) + 4 * k4 * pow(L, 2) * pow(L21, 2) + k2 * pow(L01, 2);  // 列11: 对角元
    Matrix[56] = L00 * L01 * k2 + 4 * pow(L, 2) * L10 * L11 * k3 + 4 * pow(L, 2) * L20 * L21 * k4;  // 列11: 第10行
    Matrix[57] = 6 * L * L11 * L22 * k3 - 6 * L * L12 * L21 * k4;  // 列11: 第9行
    Matrix[58] = 6 * L * L11 * L21 * k3 - 6 * L * L11 * L21 * k4;  // 列11: 第8行
    Matrix[59] = 6 * L * L11 * L20 * k3 - 6 * L * L10 * L21 * k4;  // 列11: 第7行
    Matrix[60] = 2 * pow(L, 2) * L11 * L12 * k3 - L01 * L02 * k2 + 2 * pow(L, 2) * L21 * L22 * k4;  // 列11: 第6行
    Matrix[61] = 2 * k3 * pow(L, 2) * pow(L11, 2) + 2 * k4 * pow(L, 2) * pow(L21, 2) - k2 * pow(L01, 2);  // 列11: 第5行
    Matrix[62] = 2 * pow(L, 2) * L10 * L11 * k3 - L00 * L01 * k2 + 2 * pow(L, 2) * L20 * L21 * k4;  // 列11: 第4行
    Matrix[63] = 6 * L * L12 * L21 * k4 - 6 * L * L11 * L22 * k3;  // 列11: 第3行
    Matrix[64] = 6 * L * L11 * L21 * k4 - 6 * L * L11 * L21 * k3;  // 列11: 第2行
    Matrix[65] = 6 * L * L10 * L21 * k4 - 6 * L * L11 * L20 * k3;  // 列11: 第1行

    Matrix[66] = 4 * k3 * pow(L, 2) * pow(L12, 2) + 4 * k4 * pow(L, 2) * pow(L22, 2) + k2 * pow(L02, 2);  // 列12: 对角元
    Matrix[67] = L01 * L02 * k2 + 4 * pow(L, 2) * L11 * L12 * k3 + 4 * pow(L, 2) * L21 * L22 * k4;  // 列12: 第11行
    Matrix[68] = L00 * L02 * k2 + 4 * pow(L, 2) * L10 * L12 * k3 + 4 * pow(L, 2) * L20 * L22 * k4;  // 列12: 第10行
    Matrix[69] = 6 * L * L12 * L22 * k3 - 6 * L * L12 * L22 * k4;  // 列12: 第9行
    Matrix[70] = 6 * L * L12 * L21 * k3 - 6 * L * L11 * L22 * k4;  // 列12: 第8行
    Matrix[71] = 6 * L * L12 * L20 * k3 - 6 * L * L10 * L22 * k4;  // 列12: 第7行
    Matrix[72] = 2 * k3 * pow(L, 2) * pow(L12, 2) + 2 * k4 * pow(L, 2) * pow(L22, 2) - k2 * pow(L02, 2);  // 列12: 第6行
    Matrix[73] = 2 * pow(L, 2) * L11 * L12 * k3 - L01 * L02 * k2 + 2 * pow(L, 2) * L21 * L22 * k4;  // 列12: 第5行
    Matrix[74] = 2 * pow(L, 2) * L10 * L12 * k3 - L00 * L02 * k2 + 2 * pow(L, 2) * L20 * L22 * k4;  // 列12: 第4行
    Matrix[75] = 6 * L * L12 * L22 * k4 - 6 * L * L12 * L22 * k3;  // 列12: 第3行
    Matrix[76] = 6 * L * L11 * L22 * k4 - 6 * L * L12 * L21 * k3;  // 列12: 第2行
    Matrix[77] = 6 * L * L10 * L22 * k4 - 6 * L * L12 * L20 * k3;  // 列12: 第1行
}

// Stress = axial + bending-y + bending-z + torsion shear
void CB31::ElementStress(double* stress, double* Displacement)
{
    // Get local displacements u_local = T * u_global
    double u_local[12] = {0.0};
    double T[12][12];
    CalculateTransformationMatrix(T);

    for (int i = 0; i < 12; ++i)
        for (int j = 0; j < 12; ++j)
            u_local[i] += T[i][j] * Displacement[LocationMatrix_[j] - 1];

    // Extract relevant local DOFs
    double u1x = u_local[0], u2x = u_local[6];      // axial
    double θ1y = u_local[4], θ2y = u_local[10];     // bending in Z
    double θ1z = u_local[5], θ2z = u_local[11];     // bending in Y
    double θ1x = u_local[3], θ2x = u_local[9];      // torsion (rotation about x)

    // Material and section properties
    CB31Material* mat = dynamic_cast<CB31Material*>(ElementMaterial_);
    double E = mat->E;
    double nu = mat->nu;      
    double G = mat->G;
    double Iy = mat->Iy, Iz = mat->Iz, J = mat->J;

    // Strains
    double ε_axial = (u2x - u1x) / Length_;
    double κy = (θ2z - θ1z) / Length_;
    double κz = (θ2y - θ1y) / Length_;
    double θx_prime = (θ2x - θ1x) / Length_;

    // Local coordinates of fiber point (for stress recovery)
    double y = 1.0, z = 1.0, r = sqrt(y*y + z*z);  // assume fiber at y=1, z=1
    double sigma_axial = E * (ε_axial - κy * z + κz * y);
    double τ_torsion = G * θx_prime * r;          // τ = G * θ' * r

    // Total equivalent stress (note: normal + shear not additive directly)
    *stress = sqrt(sigma_axial * sigma_axial + 3 * τ_torsion * τ_torsion); // von Mises approx.
}

void CB31::CalculateTransformationMatrix(double T[12][12])
{
    // Local x: beam axis
    double x[3] = {Direction_[0], Direction_[1], Direction_[2]};

    // Use user-provided UpVector_ as reference vector
    double up[3] = {UpVector_[0], UpVector_[1], UpVector_[2]};

    // y = up × x
    double y[3] = {
        up[1]*x[2] - up[2]*x[1],
        up[2]*x[0] - up[0]*x[2],
        up[0]*x[1] - up[1]*x[0]
    };

    // Normalize y
    double norm_y = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
    if (norm_y < 1e-10) {
        cerr << "Error: Up vector is parallel to beam axis. Cannot define local coordinate system." << endl;
        exit(-1);
    }
    for (int i = 0; i < 3; ++i) y[i] /= norm_y;

    // z = x × y
    double z[3] = {
        x[1]*y[2] - x[2]*y[1],
        x[2]*y[0] - x[0]*y[2],
        x[0]*y[1] - x[1]*y[0]
    };

    // Direction cosine matrix
    double R[3][3] = {
        {x[0], y[0], z[0]},
        {x[1], y[1], z[1]},
        {x[2], y[2], z[2]}
    };

    // Fill 12x12 T = block diagonal of R
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                T[3*i + j][3*i + k] = R[j][k];
}