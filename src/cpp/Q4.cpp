#include "Q4.h"
#include "2DMaterial.h" // 与T3共用平面材料类
#include <iostream>
#include <cmath>

using namespace std;
const int NGQP = 2;           // 高斯积分点数
const double GP[2] = { -0.577350269189626, 0.577350269189626 }; // ±1/sqrt(3)
const double W[2] = { 1.0, 1.0 };

CQ4::CQ4()
{
    NEN_ = 4;        // 4个节点
    nodes_ = new CNode * [NEN_];
    ND_ = 24;         // 4节点×6自由度
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

CQ4::~CQ4()
{
    delete[] LocationMatrix_;
}

bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet, N[4];
    Input >> N[0] >> N[1] >> N[2] >> N[3] >> MSet;

    // 材料类型
    ElementMaterial_ = dynamic_cast<C2DMaterial*>(MaterialSets) + MSet - 1 ;
   

    for (int i = 0; i < 4; i++)
        nodes_[i] = &NodeList[N[i] - 1];

    return true;
}

void CQ4::Write(COutputter& output)
{
    output << setw(9) << nodes_[0]->NodeNumber
        << setw(9) << nodes_[1]->NodeNumber
        << setw(9) << nodes_[2]->NodeNumber
        << setw(9) << nodes_[3]->NodeNumber
        << setw(12) << ElementMaterial_->nset << endl;
}

// 计算单元刚度矩阵
void CQ4::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());
    GaussIntegration(Matrix); // 通过高斯积分计算
}

// 高斯积分主函数
void CQ4::GaussIntegration(double* Matrix)
{
    C2DMaterial* material = dynamic_cast<C2DMaterial*>(ElementMaterial_);
    double E = material->E;
    double nu = material->nu;
    double t = material->t;
    int type = material->type;

    // 构造D矩阵
    double D[3][3] = { 0 };
    if (type == 0) { // 平面应力
        double factor = E / (1 - nu * nu);
        D[0][0] = factor;
        D[0][1] = factor * nu;
        D[1][1] = factor;
        D[2][2] = factor * (1 - nu) / 2;
    }
    else {         // 平面应变
        double factor = E / (1 + nu) / (1 - 2 * nu);
        D[0][0] = factor * (1 - nu);
        D[0][1] = factor * nu;
        D[1][1] = factor * (1 - nu);
        D[2][2] = factor * (1 - 2 * nu) / 2;
    }
    D[1][0] = D[0][1]; // 对称性
    double K[8][8] = { 0 };
   

    // 2x2高斯积分
    for (int i = 0; i < NGQP; ++i) {
        double xi = GP[i];
        for (int j = 0; j < NGQP; ++j) {
            double eta = GP[j];
            double N[4], dNdxi[4], dNdeta[4];
            ShapeFunction(xi, eta, N, dNdxi, dNdeta);

            double detJ;
            double dNdx[4], dNdy[4];
            Jacobian(xi, eta, detJ, dNdx, dNdy);

            // 构造B矩阵 (3x8)
            double B[3][8] = { 0 };
            ConstructBMatrix(B, dNdx, dNdy);

            double weight = W[i] * W[j] * detJ * t;

            // 计算B^T * D * B * weight
            for (int row = 0; row < 8; row++) {
                for (int col = 0; col < 8; col++) {
                    double sum = 0.0;
                    // 矩阵乘法: sum = B^T[row] * D * B[col]
                    for (int m = 0; m < 3; m++) {       // 应变分量
                        double temp = 0.0;
                        for (int n = 0; n < 3; n++) {   // D矩阵分量
                            temp += D[m][n] * B[n][col];
                        }
                        sum += B[m][row] * temp;
                    }
                    double contrib = sum * weight;
                    K[row][col] += contrib;

                    
                }
            }
        }
    }
    /*
    // 打印完整刚度矩阵
    std::cout << "Stiffness Matrix K (8x8):" << std::endl;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            std::cout << std::setw(10)
                << K[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */
    // 扩展为24x24三维格式 (z方向填充0)
    double K3D[24][24] = {0};  

    // 映射关系：每2个2D自由度 -> 3个3D自由度
    const int dofMap[8] = { 0,1,6,7,12,13,18,19 }; // 对应8个自由度的位置

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            K3D[dofMap[i]][dofMap[j]] = K[i][j];
        }

    }
   

    int index = 0; // 用于跟踪 Matrix 中的当前位置

    for (int j = 0; j < 24; j++) {
        for (int i = j; i >= 0; i--) {
            Matrix[index++] = K3D[i][j];
        }
    }
    // std::cout << "index:" << index<<std::endl;
}

      

// 构造B矩阵
void CQ4::ConstructBMatrix(double B[3][8], const double* dNdx, const double* dNdy)
{
    for (int i = 0; i < 4; ++i) {
        int posX = 2 * i;     // x自由度位置
        int posY = 2 * i + 1; // y自由度位置
        // z自由度(3*i+2)保持为0

        B[0][posX] = dNdx[i];   // ε_xx
        B[1][posY] = dNdy[i];   // ε_yy
        B[2][posX] = dNdy[i];   // γ_xy
        B[2][posY] = dNdx[i];
    }
    /*
    // 输出B矩阵
    std::cout << "B Matrix:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 8; ++j) {
            std::cout << B[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */
}

// 形函数定义
void CQ4::ShapeFunction(double xi, double eta, double* N, double* dNdxi, double* dNdeta)
{
    // 等参形函数
    if (N) {
        N[0] = 0.25 * (1 - xi) * (1 - eta);
        N[1] = 0.25 * (1 + xi) * (1 - eta);
        N[2] = 0.25 * (1 + xi) * (1 + eta);
        N[3] = 0.25 * (1 - xi) * (1 + eta);
    }
    if (dNdxi) {
        // 自然坐标导数
        dNdxi[0] = -0.25 * (1 - eta);
        dNdxi[1] = 0.25 * (1 - eta);
        dNdxi[2] = 0.25 * (1 + eta);
        dNdxi[3] = -0.25 * (1 + eta);
    }
    if (dNdeta) {
        dNdeta[0] = -0.25 * (1 - xi);
        dNdeta[1] = -0.25 * (1 + xi);
        dNdeta[2] = 0.25 * (1 + xi);
        dNdeta[3] = 0.25 * (1 - xi);
    }
}

// 雅可比矩阵计算
void CQ4::Jacobian(double xi, double eta, double& detJ, double* dNdx, double* dNdy)
{
    double dNdxi[4], dNdeta[4];
    ShapeFunction(xi, eta, nullptr, dNdxi, dNdeta);

    double J[2][2] = { 0 };
    for (int i = 0; i < 4; ++i) {
        double x = nodes_[i]->XYZ[0];
        double y = nodes_[i]->XYZ[1];
        J[0][0] += dNdxi[i] * x; // dx/dξ
        J[0][1] += dNdxi[i] * y; // dy/dξ
        J[1][0] += dNdeta[i] * x; // dx/dη
        J[1][1] += dNdeta[i] * y; // dy/dη
    }

    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    /*
    // 输出detJ
    std::cout << "detJ:" << detJ<< std::endl;
    */
    if (detJ <= 1e-10) {
        cerr << "Error: Negative Jacobian in Q4 element!" << endl;
        exit(EXIT_FAILURE);
    }

    // 计算全局导数 dNdx = dNdξ * dξ/dx + dNdη * dη/dx
    double invJ[2][2] = {
        { J[1][1] / detJ, -J[0][1] / detJ},
        {-J[1][0] / detJ,  J[0][0] / detJ}
    };

    for (int i = 0; i < 4; ++i) {
        dNdx[i] = dNdxi[i] * invJ[0][0] + dNdeta[i] * invJ[0][1];
        dNdy[i] = dNdxi[i] * invJ[1][0] + dNdeta[i] * invJ[1][1];
    }
}

// 应力计算
// 计算Q4单元4个高斯点的应力和实际坐标
void CQ4::ElementStress(double* stress, double* Displacement)
{
    // 材料参数
    C2DMaterial* material = dynamic_cast<C2DMaterial*>(ElementMaterial_);
    double E = material->E;
    double nu = material->nu;
    int type = material->type;

    // 构造D矩阵
    double D[3][3] = { 0 };
    if (type == 0) { // 平面应力
        double factor = E / (1 - nu * nu);
        D[0][0] = factor;
        D[0][1] = factor * nu;
        D[1][1] = factor;
        D[2][2] = factor * (1 - nu) / 2;
    }
    else { // 平面应变
        double factor = E / (1 + nu) / (1 - 2 * nu);
        D[0][0] = factor * (1 - nu);
        D[0][1] = factor * nu;
        D[1][1] = factor * (1 - nu);
        D[2][2] = factor * (1 - 2 * nu) / 2;
    }
    D[1][0] = D[0][1];

    // 提取二维位移
    double U[24] = { 0 };
    for (int i = 0; i < 24; ++i)
        if (LocationMatrix_[i])
            U[i] = Displacement[LocationMatrix_[i] - 1];
    double U2D[8] = { 0 };
    for (int i = 0; i < 4; ++i) {
        U2D[2 * i] = U[6 * i];
        U2D[2 * i + 1] = U[6 * i + 1];
    }

    // 2x2高斯点
    int idx = 0;
    for (int i = 0; i < 2; ++i) {
        double xi = GP[i];
        for (int j = 0; j < 2; ++j) {
            double eta = GP[j];

            // 计算B矩阵
            double dNdx[4], dNdy[4], detJ;
            Jacobian(xi, eta, detJ, dNdx, dNdy);
            double B[3][8] = { 0 };
            ConstructBMatrix(B, dNdx, dNdy);

            // 计算应变
            double strain[3] = { 0 };
            for (int m = 0; m < 3; ++m)
                for (int n = 0; n < 8; ++n)
                    strain[m] += B[m][n] * U2D[n];

            // 计算应力
            stress[idx * 3 + 0] = D[0][0] * strain[0] + D[0][1] * strain[1];
            stress[idx * 3 + 1] = D[1][0] * strain[0] + D[1][1] * strain[1];
            stress[idx * 3 + 2] = D[2][2] * strain[2];

            // 计算实际坐标
            double N[4];
            ShapeFunction(xi, eta, N, nullptr, nullptr);
            double X = 0.0, Y = 0.0;
            for (int k = 0; k < 4; ++k) {
                X += N[k] * nodes_[k]->XYZ[0];
                Y += N[k] * nodes_[k]->XYZ[1];
            }
            ++idx;
        }
    }
}

// 计算4个高斯点的实际坐标
void CQ4::ElementGauss(double* gp_coords)
{
    // 2x2高斯点
   
    int idx = 0;
    for (int i = 0; i < 2; ++i) {
        double xi = GP[i];
        for (int j = 0; j < 2; ++j) {
            double eta = GP[j];
            double N[4];
            ShapeFunction(xi, eta, N, nullptr, nullptr);
            double X = 0.0, Y = 0.0;
            for (int k = 0; k < 4; ++k) {
                X += N[k] * nodes_[k]->XYZ[0];
                Y += N[k] * nodes_[k]->XYZ[1];
            }
            gp_coords[idx * 2 + 0] = X;
            gp_coords[idx * 2 + 1] = Y;
            ++idx;
        }
    }
}


