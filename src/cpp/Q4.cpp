#include "Q4.h"
#include "2DMaterial.h" // ��T3����ƽ�������
#include <iostream>
#include <cmath>

using namespace std;
const int NGQP = 2;           // ��˹���ֵ���
const double GP[2] = { -0.577350269189626, 0.577350269189626 }; // ��1/sqrt(3)
const double W[2] = { 1.0, 1.0 };

CQ4::CQ4()
{
    NEN_ = 4;        // 4���ڵ�
    nodes_ = new CNode * [NEN_];
    ND_ = 24;         // 4�ڵ��6���ɶ�
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

    // ��������
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

// ���㵥Ԫ�նȾ���
void CQ4::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());
    GaussIntegration(Matrix); // ͨ����˹���ּ���
}

// ��˹����������
void CQ4::GaussIntegration(double* Matrix)
{
    C2DMaterial* material = dynamic_cast<C2DMaterial*>(ElementMaterial_);
    double E = material->E;
    double nu = material->nu;
    double t = material->t;
    int type = material->type;

    // ����D����
    double D[3][3] = { 0 };
    if (type == 0) { // ƽ��Ӧ��
        double factor = E / (1 - nu * nu);
        D[0][0] = factor;
        D[0][1] = factor * nu;
        D[1][1] = factor;
        D[2][2] = factor * (1 - nu) / 2;
    }
    else {         // ƽ��Ӧ��
        double factor = E / (1 + nu) / (1 - 2 * nu);
        D[0][0] = factor * (1 - nu);
        D[0][1] = factor * nu;
        D[1][1] = factor * (1 - nu);
        D[2][2] = factor * (1 - 2 * nu) / 2;
    }
    D[1][0] = D[0][1]; // �Գ���
    double K[8][8] = { 0 };
   

    // 2x2��˹����
    for (int i = 0; i < NGQP; ++i) {
        double xi = GP[i];
        for (int j = 0; j < NGQP; ++j) {
            double eta = GP[j];
            double N[4], dNdxi[4], dNdeta[4];
            ShapeFunction(xi, eta, N, dNdxi, dNdeta);

            double detJ;
            double dNdx[4], dNdy[4];
            Jacobian(xi, eta, detJ, dNdx, dNdy);

            // ����B���� (3x8)
            double B[3][8] = { 0 };
            ConstructBMatrix(B, dNdx, dNdy);

            double weight = W[i] * W[j] * detJ * t;

            // ����B^T * D * B * weight
            for (int row = 0; row < 8; row++) {
                for (int col = 0; col < 8; col++) {
                    double sum = 0.0;
                    // ����˷�: sum = B^T[row] * D * B[col]
                    for (int m = 0; m < 3; m++) {       // Ӧ�����
                        double temp = 0.0;
                        for (int n = 0; n < 3; n++) {   // D�������
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
    // ��ӡ�����նȾ���
    std::cout << "Stiffness Matrix K (8x8):" << std::endl;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            std::cout << std::setw(10)
                << K[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */
    // ��չΪ24x24��ά��ʽ (z�������0)
    double K3D[24][24] = {0};  

    // ӳ���ϵ��ÿ2��2D���ɶ� -> 3��3D���ɶ�
    const int dofMap[8] = { 0,1,6,7,12,13,18,19 }; // ��Ӧ8�����ɶȵ�λ��

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            K3D[dofMap[i]][dofMap[j]] = K[i][j];
        }

    }
   

    int index = 0; // ���ڸ��� Matrix �еĵ�ǰλ��

    for (int j = 0; j < 24; j++) {
        for (int i = j; i >= 0; i--) {
            Matrix[index++] = K3D[i][j];
        }
    }
    // std::cout << "index:" << index<<std::endl;
}

      

// ����B����
void CQ4::ConstructBMatrix(double B[3][8], const double* dNdx, const double* dNdy)
{
    for (int i = 0; i < 4; ++i) {
        int posX = 2 * i;     // x���ɶ�λ��
        int posY = 2 * i + 1; // y���ɶ�λ��
        // z���ɶ�(3*i+2)����Ϊ0

        B[0][posX] = dNdx[i];   // ��_xx
        B[1][posY] = dNdy[i];   // ��_yy
        B[2][posX] = dNdy[i];   // ��_xy
        B[2][posY] = dNdx[i];
    }
    /*
    // ���B����
    std::cout << "B Matrix:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 8; ++j) {
            std::cout << B[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */
}

// �κ�������
void CQ4::ShapeFunction(double xi, double eta, double* N, double* dNdxi, double* dNdeta)
{
    // �Ȳ��κ���
    if (N) {
        N[0] = 0.25 * (1 - xi) * (1 - eta);
        N[1] = 0.25 * (1 + xi) * (1 - eta);
        N[2] = 0.25 * (1 + xi) * (1 + eta);
        N[3] = 0.25 * (1 - xi) * (1 + eta);
    }
    if (dNdxi) {
        // ��Ȼ���굼��
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

// �ſɱȾ������
void CQ4::Jacobian(double xi, double eta, double& detJ, double* dNdx, double* dNdy)
{
    double dNdxi[4], dNdeta[4];
    ShapeFunction(xi, eta, nullptr, dNdxi, dNdeta);

    double J[2][2] = { 0 };
    for (int i = 0; i < 4; ++i) {
        double x = nodes_[i]->XYZ[0];
        double y = nodes_[i]->XYZ[1];
        J[0][0] += dNdxi[i] * x; // dx/d��
        J[0][1] += dNdxi[i] * y; // dy/d��
        J[1][0] += dNdeta[i] * x; // dx/d��
        J[1][1] += dNdeta[i] * y; // dy/d��
    }

    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    /*
    // ���detJ
    std::cout << "detJ:" << detJ<< std::endl;
    */
    if (detJ <= 1e-10) {
        cerr << "Error: Negative Jacobian in Q4 element!" << endl;
        exit(EXIT_FAILURE);
    }

    // ����ȫ�ֵ��� dNdx = dNd�� * d��/dx + dNd�� * d��/dx
    double invJ[2][2] = {
        { J[1][1] / detJ, -J[0][1] / detJ},
        {-J[1][0] / detJ,  J[0][0] / detJ}
    };

    for (int i = 0; i < 4; ++i) {
        dNdx[i] = dNdxi[i] * invJ[0][0] + dNdeta[i] * invJ[0][1];
        dNdy[i] = dNdxi[i] * invJ[1][0] + dNdeta[i] * invJ[1][1];
    }
}

// Ӧ������
// ����Q4��Ԫ4����˹���Ӧ����ʵ������
void CQ4::ElementStress(double* stress, double* Displacement)
{
    // ���ϲ���
    C2DMaterial* material = dynamic_cast<C2DMaterial*>(ElementMaterial_);
    double E = material->E;
    double nu = material->nu;
    int type = material->type;

    // ����D����
    double D[3][3] = { 0 };
    if (type == 0) { // ƽ��Ӧ��
        double factor = E / (1 - nu * nu);
        D[0][0] = factor;
        D[0][1] = factor * nu;
        D[1][1] = factor;
        D[2][2] = factor * (1 - nu) / 2;
    }
    else { // ƽ��Ӧ��
        double factor = E / (1 + nu) / (1 - 2 * nu);
        D[0][0] = factor * (1 - nu);
        D[0][1] = factor * nu;
        D[1][1] = factor * (1 - nu);
        D[2][2] = factor * (1 - 2 * nu) / 2;
    }
    D[1][0] = D[0][1];

    // ��ȡ��άλ��
    double U[24] = { 0 };
    for (int i = 0; i < 24; ++i)
        if (LocationMatrix_[i])
            U[i] = Displacement[LocationMatrix_[i] - 1];
    double U2D[8] = { 0 };
    for (int i = 0; i < 4; ++i) {
        U2D[2 * i] = U[6 * i];
        U2D[2 * i + 1] = U[6 * i + 1];
    }

    // 2x2��˹��
    int idx = 0;
    for (int i = 0; i < 2; ++i) {
        double xi = GP[i];
        for (int j = 0; j < 2; ++j) {
            double eta = GP[j];

            // ����B����
            double dNdx[4], dNdy[4], detJ;
            Jacobian(xi, eta, detJ, dNdx, dNdy);
            double B[3][8] = { 0 };
            ConstructBMatrix(B, dNdx, dNdy);

            // ����Ӧ��
            double strain[3] = { 0 };
            for (int m = 0; m < 3; ++m)
                for (int n = 0; n < 8; ++n)
                    strain[m] += B[m][n] * U2D[n];

            // ����Ӧ��
            stress[idx * 3 + 0] = D[0][0] * strain[0] + D[0][1] * strain[1];
            stress[idx * 3 + 1] = D[1][0] * strain[0] + D[1][1] * strain[1];
            stress[idx * 3 + 2] = D[2][2] * strain[2];

            // ����ʵ������
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

// ����4����˹���ʵ������
void CQ4::ElementGauss(double* gp_coords)
{
    // 2x2��˹��
   
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


