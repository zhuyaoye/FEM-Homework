// CKirchhoffPlate.cpp
#include "CKirchhoffPlate.h"
#include <iostream>
#include <iomanip>
#include <cmath>


using namespace std;

// Constructor
CKirchhoffPlate::CKirchhoffPlate()
{
    NEN_ = 4;    // 4 nodes per plate element
    nodes_ = new CNode*[NEN_];
    ND_ = 24;    // 3 DOF per node
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}

// Destructor
CKirchhoffPlate::~CKirchhoffPlate()
{
    delete[] nodes_;
    delete[] LocationMatrix_;
}

// Read element data from input
bool CKirchhoffPlate::Read(ifstream& input, CMaterial* materialSets, CNode* nodeList)
{
    unsigned int nodeNumbers[4], mset;
    for (int i = 0; i < 4; ++i)
        input >> nodeNumbers[i];
    input >> mset;

    for (int i = 0; i < 4; ++i)
        nodes_[i] = &nodeList[nodeNumbers[i] - 1];

    ElementMaterial_ = dynamic_cast<CPlateMaterial*>(materialSets + mset - 1);
    return true;
}

// Write element data to output
void CKirchhoffPlate::Write(COutputter& output)
{
    for (int i = 0; i < 4; ++i)
        output << setw(9) << nodes_[i]->NodeNumber;
    output << setw(12) << ElementMaterial_->nset << endl;
}



// Compute stiffness matrix using 2x2 Gauss quadrature
void CKirchhoffPlate::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(12, 12);      // Kirchhoff 板单元刚度矩阵
    Eigen::MatrixXd ke_full  = Eigen::MatrixXd::Zero(24, 24);      // 扩展后用于系统装配的刚度矩阵
    Eigen::Matrix3d D = ElementMaterial_->GetFlexuralMatrix();

    double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
    
    double a = sqrt(DX[0] *  DX[0] + DX[1] * DX[1] + DX[2] * DX[2])/2;

    for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[3]->XYZ[i] - nodes_[0]->XYZ[i];
    
    double b = sqrt(DX[0] *  DX[0] + DX[1] * DX[1] + DX[2] * DX[2])/2;


    double t = ElementMaterial_->GetThickness();
    double J = a * b ;  // constant for rectangular element

    const double gp[2] = { -sqrt(1.0 / 3.0), sqrt(1.0 / 3.0) };
    const double w[2] = { 1.0, 1.0 };

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            double xi = gp[i];
            double eta = gp[j];
            Eigen::Matrix<double, 3, 12> B = ComputeBMatrix(xi, eta, a, b);
            ke_local += B.transpose() * D * B * 1 * J * w[i] * w[j];
            //std::cout << "K matrix (12x12):\n" << ke_local << "\n\n";
            //std::cout << "B matrix (12x12):\n" << B << "\n\n";
            //std::cout << "D matrix (12x12):\n" << D << "\n\n";
        }
    }
     // DOF 映射表：原来 12 DOF → 扩展系统 24 DOF
    int dof_map[12] = {
        2, 3, 4,     // node 0: w, θx, θy → global DOF 2,3,4
        8, 9,10,     // node 1: ...
       14,15,16,     // node 2
       20,21,22      // node 3
    };

    // 将 ke_local 嵌入到 ke_full 中
    for (int i = 0; i < 12; ++i)
    {
        for (int j = 0; j < 12; ++j)
        {
            int row = dof_map[i];
            int col = dof_map[j];
            ke_full(row, col) += ke_local(i, j);  // 注意 += 可兼容重复贡献
        }
    }
    //std::cout << "K matrix (12x12):\n" << ke_full << "\n\n";
    // 打包 ke_full (24x24) 上三角到 Matrix[]
    int index = 0;
    const int n = 24;
    for (int k = 0; k < n; ++k)
    {
        Matrix[index++] = ke_full(k, k);
        for (int i = k - 1; i >= 0; --i)
        {
            Matrix[index++] = ke_full(i, k);
        }
    }
}

// Compute B-matrix at a given Gauss point
Eigen::Matrix<double, 3, 12> CKirchhoffPlate::ComputeBMatrix(double xi, double eta, double a, double b)
{
    Eigen::Matrix<double, 3, 12> B = Eigen::Matrix<double, 3, 12>::Zero();

    double ab = 4.0 * a * b;
    
    // Fill B matrix using same expressions as your original code
    B(0,0) = 3.0 * b/a * xi * (1 - eta) / ab;
    B(0,2) = b * (1 - 3.0 * xi) * (1 - eta) / ab;
    B(0,3) = -3.0 * b/a * xi * (1 - eta) / ab;
    B(0,5) = -b * (1 + 3.0 * xi) * (1 - eta) / ab;
    B(0,6) = -3.0 * b/a * xi * (1 + eta) / ab;
    B(0,8) = -b * (1 + 3.0 * xi) * (1 + eta) / ab;
    B(0,9) = 3.0 * b/a * xi * (1 + eta) / ab;
    B(0,11)= b  * (1 - 3.0 * xi) * (1 + eta) / ab;

    B(1,0) = 3.0 * a/b * eta * (1 - xi) / ab;
    B(1,1) = -a * (1 - 3.0 * eta) * (1 - xi) / ab;
    B(1,3) = 3.0 * a/b * eta * (1 + xi) / ab;
    B(1,4) = -a  * (1 - 3.0 * eta) * (1 + xi) / ab;
    B(1,6) = -3.0 * a/b * eta * (1 + xi) / ab;
    B(1,7) = a * (1 + 3.0 * eta) * (1 + xi) / ab;
    B(1,9)= -3.0 * a/b * eta * (1 - xi) / ab;
    B(1,10)= a * (1 + 3.0 * eta) * (1 - xi) / ab;

    
    B(2, 0) = (4.0 - 3.0 * xi * xi - 3.0 * eta * eta) / ab;
    B(2, 1) = -b * (3.0 * eta * eta - 2.0 * eta- 1.0) / ab;
    B(2, 2) = -a  * (1.0 + 2.0  * xi - 3.0 * xi * xi) / ab;
    B(2, 3) = -(4.0 - 3.0 * xi * xi - 3.0 * eta * eta) / ab;
    B(2, 4) = b * (3.0 * eta * eta - 2.0 * eta- 1.0) / ab;
    B(2, 5) = -a  * (1.0 - 2.0  * xi - 3.0 * xi * xi) / ab;
    B(2, 6) = (4.0 - 3.0 * xi * xi - 3.0 * eta * eta) / ab;
    B(2, 7) = b * (3.0 * eta * eta + 2.0 * eta- 1.0) / ab;
    B(2, 8) = a  * (1.0 - 2.0  * xi - 3.0 * xi * xi) / ab;
    B(2, 9) = -(4.0 - 3.0 * xi * xi - 3.0 * eta * eta) / ab;
    B(2, 10) = -b * (3.0 * eta * eta + 2.0 * eta- 1.0) / ab;
    B(2, 11) = a  * (1.0 + 2.0  * xi - 3.0 * xi * xi) / ab;
    

    return B;
}

// Compute equivalent load vector
void CKirchhoffPlate::ElementLoad(double* Load, double q)
{
    Eigen::VectorXd fe = Eigen::VectorXd::Zero(12);
    
    double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
    
    double a = sqrt(DX[0] *  DX[0] + DX[1] * DX[1] + DX[2] * DX[2])/2;

    for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[3]->XYZ[i] - nodes_[0]->XYZ[i];
    
    double b = sqrt(DX[0] *  DX[0] + DX[1] * DX[1] + DX[2] * DX[2])/2;

    Eigen::Vector3d f1 = (q * a * b / 3.0) * Eigen::Vector3d(3, b, -a);
    Eigen::Vector3d f2 = (q * a * b / 3.0) * Eigen::Vector3d(3,  -b, a);
    Eigen::Vector3d f3 = (q * a * b / 3.0) * Eigen::Vector3d(3,  b,  a);
    Eigen::Vector3d f4 = (q * a * b / 3.0) * Eigen::Vector3d(3, -b,  a);

    fe.segment<3>(0)  = f1;
    fe.segment<3>(3)  = f2;
    fe.segment<3>(6)  = f3;
    fe.segment<3>(9)  = f4;

    for (int i = 0; i < 12; ++i)
        Load[i] = fe(i);
}

void CKirchhoffPlate::ElementStress(double* stress, double* Displacement)
{
    Eigen::Matrix3d D = ElementMaterial_->GetFlexuralMatrix();   //用于计算M（中面）
    //Eigen::Matrix3d D = ElementMaterial_->GetstressMatrix();  //可用于计算应力
    // 获取矩形单元的长宽
    double DX[3];
    for (int i = 0; i < 3; ++i)
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
    double a = sqrt(DX[0]*DX[0] + DX[1]*DX[1] + DX[2]*DX[2])/2;

    for (int i = 0; i < 3; ++i)
        DX[i] = nodes_[3]->XYZ[i] - nodes_[0]->XYZ[i];
    double b = sqrt(DX[0]*DX[0] + DX[1]*DX[1] + DX[2]*DX[2])/2;

    const double gp[2] = { -sqrt(1.0 / 3.0), sqrt(1.0 / 3.0) };
    //const double gp[2] = { -1, 1 };
    
    // 预留空间存储 4 个 Gauss 点的中面上广义内力 (每点 3 个分量)

    int ind = 0;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            double xi = gp[i];
            double eta = gp[j];

            // 构造 B 矩阵
            Eigen::Matrix<double, 3, 12> B = ComputeBMatrix(xi, eta, a, b);
            
            // 将位移拷贝成 Eigen 向量
            Eigen::VectorXd de(12);
            for (int node = 0; node < 4; ++node)
            {
                int base = node * 6;  // 每节点 6 个自由度
                int lm_w  = LocationMatrix_[base + 2];  // w
                int lm_tx = LocationMatrix_[base + 3];  // θx
                int lm_ty = LocationMatrix_[base + 4];  // θy

                de(node * 3 + 0) = (lm_w  > 0) ? Displacement[lm_w - 1]  : 0.0;
                de(node * 3 + 1) = (lm_tx > 0) ? Displacement[lm_tx - 1] : 0.0;
                de(node * 3 + 2) = (lm_ty > 0) ? Displacement[lm_ty - 1] : 0.0;
            }


            // 计算应变 (曲率)
            Eigen::Vector3d strain = B * de;

            double z_up = ElementMaterial_->GetThickness() / 2.0;  // 上表面 z = +t/2
            Eigen::Vector3d sigma =  -D * strain;
            //Eigen::Vector3d sigma =  -D * z_up * strain;  //求应力，记得修改D（和求M的D差常数）
            // 存储结果
            for (int s = 0; s < 3; ++s)
                stress[3 * ind + s] = sigma(s);
            ++ind;
        }
    }
}