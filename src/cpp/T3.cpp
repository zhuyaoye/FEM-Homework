#include "T3.h"
#include "2DMaterial.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

CT3::CT3()
{
    NEN_ = 3;       // 3 nodes
    nodes_ = new CNode * [NEN_];
    ND_ = 18;        // 6 DOF per node * 3
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}

CT3::~CT3()
{
}

bool CT3::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet, N[3];
    Input >> N[0] >> N[1] >> N[2] >> MSet;

    for (int i = 0; i < 3; i++)
        nodes_[i] = &NodeList[N[i] - 1];

    ElementMaterial_ = dynamic_cast<C2DMaterial*>(MaterialSets) + MSet - 1;
    return true;
}

void CT3::Write(COutputter& output)
{
    output << setw(9) << nodes_[0]->NodeNumber
        << setw(9) << nodes_[1]->NodeNumber
        << setw(9) << nodes_[2]->NodeNumber
        << setw(12) << ElementMaterial_->nset << endl;
}

// element stiffness matrux
void CT3::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    C2DMaterial* material_ = dynamic_cast<C2DMaterial*>(ElementMaterial_);	// Pointer to material of the element
    
    // Material properties
    double E = material_->E;
    double nu = material_->nu;
    double t = material_->t;
    int type = material_->type;

    // ���㵥Ԫ���
    double A = CalculateArea();

    // �����κ�������
    double dNdx[3], dNdy[3];
    ShapeFunction(nullptr, dNdx, dNdy);

    // ����B����
    double B[3][6] = { 0 };
    for (int i = 0; i < 3; i++) {
        B[0][2 * i] = dNdx[i];
        B[1][2 * i + 1] = dNdy[i];
        B[2][2 * i] = dNdy[i];
        B[2][2 * i + 1] = dNdx[i];
    }
    /*
    // ��ӡB����
    std::cout << "Matrix B (3x6):" << std::endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 6; ++j) {
            std::cout << std::setw(10)
                << B[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */
    // ����D����
    double D[3][3] = { 0 };
    if (type == 0) { // ƽ��Ӧ��
        double factor = E / (1 - nu * nu);
        D[0][0] = factor;
        D[0][1] = factor * nu;
        D[1][0] = factor * nu;
        D[1][1] = factor;
        D[2][2] = factor * (1 - nu) / 2;
    }
    else {       // ƽ��Ӧ��
        double factor = E / (1 + nu) / (1 - 2 * nu);
        D[0][0] = factor * (1 - nu);
        D[0][1] = factor * nu;
        D[1][0] = factor * nu;
        D[1][1] = factor * (1 - nu);
        D[2][2] = factor * (1 - 2 * nu) / 2;
    }

	double K[6][6] = { 0 }; // �նȾ���
    // ����նȾ��� K = B^T*D*B*t*A
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 3; k++)
                K[i][j]+= ( B[k][i] * D[k][0] * B[0][j] +
                B[k][i] * D[k][1] * B[1][j] +
                B[k][i] * D[k][2] * B[2][j]) * t * A;
        }
    }

    /*
   // ��ӡ�����նȾ���
   std::cout << "Stiffness Matrix K (6x6):" << std::endl;
   for (int i = 0; i < 6; ++i) {
       for (int j = 0; j < 6; ++j) {
           std::cout << std::setw(10)
               << K[i][j] << " ";
       }
       std::cout << std::endl;
   }
   */
   // ��չΪ18*18��ά��ʽ (z�������0)
    double K3D[18][18] = { 0 };

    // ӳ���ϵ��ÿ2��2D���ɶ� -> 6��3D���ɶ�
    const int dofMap[6] = { 0,1,6,7,12,13 }; // ��Ӧ8�����ɶȵ�λ��

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            K3D[dofMap[i]][dofMap[j]] = K[i][j];
        }

    }
    int index = 0; // ���ڸ��� Matrix �еĵ�ǰλ��

    for (int j = 0; j < 18; j++) {
        for (int i = j; i >= 0; i--) {
            Matrix[index++] = K3D[i][j];
        }
    }
    // std::cout << "index:" << index<<std::endl;
}

// ���㵥ԪӦ��
void CT3::ElementStress(double* stress, double* Displacement)
{
    // ��ȡ���ϲ���
    C2DMaterial* material_ = dynamic_cast<C2DMaterial*>(ElementMaterial_);
    double E = material_->E;
    double nu = material_->nu;
    int type = material_->type;

    // �����κ�������
    double dNdx[3], dNdy[3];
    ShapeFunction(nullptr, dNdx, dNdy);

    // ����B����
    double B[3][6] = { 0 };
    for (int i = 0; i < 3; i++) {
        B[0][2 * i] = dNdx[i];
        B[1][2 * i + 1] = dNdy[i];
        B[2][2 * i] = dNdy[i];
        B[2][2 * i + 1] = dNdx[i];
    }

    // ��ȡ��Ԫλ��
    double U[18] = {0};
    for (int i = 0; i < 18; i++)
        if (LocationMatrix_[i])
        U[i] = Displacement[LocationMatrix_[i] - 1];

    // ��ȡ��άλ�ƣ�����z�������ɶȣ�
    double U2D[6] = { 0 };  // ��ά������Ҫ8�����ɶ� (4�ڵ� x 2���ɶ�)
    for (int i = 0; i < 3; ++i) {
        // ÿ���ڵ�ȡǰ�������ɶ� (x,y), ����������(z)

        U2D[2 * i] = U[6* i];    // x����λ��
        U2D[2 * i + 1] = U[6 * i + 1];  // y����λ��
    }
    // ����Ӧ�� �� = B*U
    double strain[3] = { 0 };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 6; j++)
            strain[i] += B[i][j] * U2D[j];

    // ����Ӧ�� �� = D*��
    if (type == 0) { // ƽ��Ӧ��
        stress[0] = E / (1 - nu * nu) * (strain[0] + nu * strain[1]);
        stress[1] = E / (1 - nu * nu) * (nu * strain[0] + strain[1]);
        stress[2] = E / (2 * (1 + nu)) * strain[2];
    }
    else {       // ƽ��Ӧ��
        stress[0] = E / (1 + nu) / (1 - 2 * nu) * ((1 - nu) * strain[0] + nu * strain[1]);
        stress[1] = E / (1 + nu) / (1 - 2 * nu) * (nu * strain[0] + (1 - nu) * strain[1]);
        stress[2] = E / (2 * (1 + nu)) * strain[2];
    }
}

// ���㵥Ԫ���
double CT3::CalculateArea()
{
    double x1 = nodes_[0]->XYZ[0], y1 = nodes_[0]->XYZ[1];
    double x2 = nodes_[1]->XYZ[0], y2 = nodes_[1]->XYZ[1];
    double x3 = nodes_[2]->XYZ[0], y3 = nodes_[2]->XYZ[1];

    return 0.5 * fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
}

// �����κ�����������Ӧ�䵥Ԫ��
void CT3::ShapeFunction(double* N, double* dNdx, double* dNdy)
{
    double x1 = nodes_[0]->XYZ[0], y1 = nodes_[0]->XYZ[1];
    double x2 = nodes_[1]->XYZ[0], y2 = nodes_[1]->XYZ[1];
    double x3 = nodes_[2]->XYZ[0], y3 = nodes_[2]->XYZ[1];

    double detJ = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

    dNdx[0] = (y2 - y3) / detJ;
    dNdx[1] = (y3 - y1) / detJ;
    dNdx[2] = (y1 - y2) / detJ;

    dNdy[0] = (x3 - x2) / detJ;
    dNdy[1] = (x1 - x3) / detJ;
    dNdy[2] = (x2 - x1) / detJ;
}