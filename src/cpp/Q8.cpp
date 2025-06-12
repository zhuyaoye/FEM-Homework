#include "Q8.h"
#include "2DMaterial.h"
#include <iostream>
#include <cmath>


const double GP_Reduced[2] = { -sqrt(1.0 / 3.0), sqrt(1.0 / 3.0) };
const double W_Reduced[2] = { 1.0, 1.0 };
const int NGQP = 2;

CQ8::CQ8() 
{
    NEN_ = 8;        // 8���ڵ�
    nodes_ = new CNode * [NEN_];
    ND_ = 48;         // 8�ڵ��6���ɶ�
    LocationMatrix_ = new unsigned int[ND_];
}

CQ8::~CQ8()
{
    delete[] LocationMatrix_;
}

bool CQ8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet, N[8];
    Input >> N[0] >> N[1] >> N[2] >> N[3] >> N[4] >> N[5] >> N[6] >> N[7] >> MSet;

    ElementMaterial_ = dynamic_cast<C2DMaterial*>(MaterialSets) + MSet - 1;
   
    for (int i = 0; i < 8; i++)
        nodes_[i] = &NodeList[N[i] - 1];

    return true;
}

void CQ8::Write(COutputter& output)
{
    output << setw(5) << nodes_[0]->NodeNumber
        << setw(5) << nodes_[1]->NodeNumber
        << setw(5) << nodes_[2]->NodeNumber
        << setw(5) << nodes_[3]->NodeNumber
        << setw(5) << nodes_[4]->NodeNumber
        << setw(5) << nodes_[5]->NodeNumber
        << setw(5) << nodes_[6]->NodeNumber
        << setw(5) << nodes_[7]->NodeNumber
        << setw(8) << ElementMaterial_->nset << endl;
}

// �նȾ���ͨ����˹���ּ���
void CQ8::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());
    GaussIntegration(Matrix);
}

// ��˹���ֺ���
void CQ8::GaussIntegration(double* Matrix)
{

    C2DMaterial* material_ = dynamic_cast<C2DMaterial*>(ElementMaterial_);
    double E = material_->E;
    double nu = material_->nu;
    double t = material_->t;
    int type = material_->type;

    // ȷ�����ֵ���Ȩ��
    const double* GP_X, * GP_Y, * W_X, * W_Y;


    GP_X = GP_Reduced; GP_Y = GP_Reduced;
    W_X = W_Reduced; W_Y = W_Reduced;


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
    double K[16][16] = { 0 };
    // ��˹����ѭ��
    for (int i = 0; i < NGQP; ++i) {
        double xi = GP_X[i];
        for (int j = 0; j < NGQP; ++j) {
            double eta = GP_Y[j];

            double N[8], dNdxi[8], dNdeta[8];
            ShapeFunction(xi, eta, N, dNdxi, dNdeta);

            double dNdx[8], dNdy[8], detJ;
            Jacobian(xi, eta, detJ, dNdx, dNdy);

            double B[3][16] = { 0 };  // ��ʼ��Ϊ0
            for (int k = 0; k < 8; k++) {
                int col = 2 * k;
                B[0][col] = dNdx[k];   // ��_xx
                B[1][col + 1] = dNdy[k]; // ��_yy
                B[2][col] = dNdy[k];   // ��_xy
                B[2][col + 1] = dNdx[k];
             
            }

            double weight = W_X[i] * W_Y[j] * detJ * t;
            // ����B^T * D * B * weight
            for (int row = 0; row < 16; row++) {
                for (int col = 0; col < 16; col++) {
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

    // ��չΪ12x12��ά��ʽ (z�������0)
    double K3D[48][48] = { 0 };

    // ӳ���ϵ��ÿ2��2D���ɶ� -> 6��3D���ɶ�
    const int dofMap[16] = { 0,1, 6,7, 12,13, 18,19, 24,25, 30,31, 36,37, 42,43 }; // ��Ӧ8�����ɶȵ�λ��

    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            K3D[dofMap[i]][dofMap[j]] = K[i][j];
        }

    }


    int index = 0; // ���ڸ��� Matrix �еĵ�ǰλ��

    for (int j = 0; j < 48; j++) {
        for (int i = j; i >= 0; i--) {
            Matrix[index++] = K3D[i][j];
        }
    }
}

// �˽ڵ��κ�����Serendipity��
void CQ8::ShapeFunction(double xi, double eta, double* N, double* dNdxi, double* dNdeta)
{
    if (N) {
        N[0] = 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1);
        N[1] = 0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1);
        N[2] = 0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1);
        N[3] = 0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1);

        // ���нڵ�
        N[4] = 0.5 * (1 - xi * xi) * (1 - eta);
        N[5] = 0.5 * (1 + xi) * (1 - eta * eta);
        N[6] = 0.5 * (1 - xi * xi) * (1 + eta);
        N[7] = 0.5 * (1 - xi) * (1 - eta * eta);
    }
    // ��Ȼ���굼��
    if (dNdxi) {
        dNdxi[0] = -0.25 * (1 - eta) * (-xi - eta - 1) - 0.25 * (1 - xi) * (1 - eta);
        dNdxi[1] = 0.25 * (1 - eta) * (xi - eta - 1) + 0.25 * (1 + xi) * (1 - eta);
        dNdxi[2] = 0.25 * (1 + xi) * (1 + eta) + 0.25 * (1 + eta) * (xi + eta - 1);
        dNdxi[3] = -0.25 * (1 + eta) * (-xi + eta - 1) - 0.25 * (1 - xi) * (1 + eta);
        dNdxi[4] = -xi * (1 - eta);
        dNdxi[5] = 0.5 * (1 - eta * eta);
        dNdxi[6] = -xi * (1 + eta);
        dNdxi[7] = -0.5 * (1 - eta * eta);
    }
    if (dNdeta) {
        dNdeta[0] = -0.25 * (1 - xi) * (1 - eta) - 0.25 * (1 - xi) * (-xi - eta - 1);
        dNdeta[1] = -0.25 * (1 + xi) * (1 - eta) - 0.25 * (1 + xi) * (xi - eta - 1);
        dNdeta[2] = 0.25 * (1 + xi) * (xi + eta - 1) + 0.25 * (1 + eta) * (1 + xi);
        dNdeta[3] = 0.25 * (1 - xi) * (1 + eta) + 0.25 * (1 - xi) * (-xi + eta - 1);
        dNdeta[4] = -0.5 * (1 - xi * xi);
        dNdeta[5] = -eta * (1 + xi);
        dNdeta[6] = 0.5 * (1 - xi * xi);
        dNdeta[7] = -eta * (1 - xi);
    }
}

// �ſɱȾ�����Q4���ƣ����ڵ�����ͬ��
void CQ8::Jacobian(double xi, double eta, double& detJ, double* dNdx, double* dNdy)
{
    double dNdxi[8], dNdeta[8];
    ShapeFunction(xi, eta, nullptr, dNdxi, dNdeta);

    double J[2][2] = { 0 };
    for (int i = 0; i < 8; ++i) {
        double x = nodes_[i]->XYZ[0];
        double y = nodes_[i]->XYZ[1];
        J[0][0] += dNdxi[i] * x;
        J[0][1] += dNdxi[i] * y;
        J[1][0] += dNdeta[i] * x;
        J[1][1] += dNdeta[i] * y;
    }

    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    
    // ���detJ
    //std::cout << "detJ:" << detJ<< std::endl;
    
    if (detJ <= 1e-10) {
        cerr << "Error: Negative Jacobian in Q8 element!" << endl;
        exit(EXIT_FAILURE);
    }

    // ��������
    double invJ[2][2] = {
        { J[1][1] / detJ, -J[0][1] / detJ},
        {-J[1][0] / detJ,  J[0][0] / detJ}
    };

    for (int i = 0; i < 8; ++i) {
        dNdx[i] = dNdxi[i] * invJ[0][0] + dNdeta[i] * invJ[0][1];
        dNdy[i] = dNdxi[i] * invJ[1][0] + dNdeta[i] * invJ[1][1];
    }
}



// Ӧ������
// ����4����˹���Ӧ��
void CQ8::ElementStress(double* stress, double* Displacement)
{
    // 2x2��˹��
    int idx = 0;
    for (int i = 0; i < 2; ++i) {
        double xi = GP_Reduced[i];
        for (int j = 0; j < 2; ++j) {
            double eta = GP_Reduced[j];

            // 1. ����Ӧ��
            double dNdx[8], dNdy[8], detJ;
            Jacobian(xi, eta, detJ, dNdx, dNdy);

            double B[3][16] = { 0 };
            for (int k = 0; k < 8; k++) {
                int col = 2 * k;
                B[0][col] = dNdx[k];
                B[1][col + 1] = dNdy[k];
                B[2][col] = dNdy[k];
                B[2][col + 1] = dNdx[k];
            }

            double strain[3] = { 0 };
            double U[48] = { 0 };
            for (int m = 0; m < 48; ++m)
                if (LocationMatrix_[m])
                    U[m] = Displacement[LocationMatrix_[m] - 1];

            double U2D[16] = { 0 };
            for (int m = 0; m < 8; ++m) {
                U2D[2 * m] = U[6 * m];
                U2D[2 * m + 1] = U[6 * m + 1];
            }
            for (int m = 0; m < 3; ++m)
                for (int n = 0; n < 16; ++n)
                    strain[m] += B[m][n] * U2D[n];

            // 2. ����Ӧ��
            C2DMaterial* material_ = dynamic_cast<C2DMaterial*>(ElementMaterial_);
            double E = material_->E, nu = material_->nu;
            int type = material_->type;

            if (type == 0) { // ƽ��Ӧ��
                stress[idx * 3 + 0] = E / (1 - nu * nu) * (strain[0] + nu * strain[1]);
                stress[idx * 3 + 1] = E / (1 - nu * nu) * (strain[1] + nu * strain[0]);
                stress[idx * 3 + 2] = E / (2 * (1 + nu)) * strain[2];
            }
            else {        // ƽ��Ӧ��
                double factor = E / ((1 + nu) * (1 - 2 * nu));
                stress[idx * 3 + 0] = factor * ((1 - nu) * strain[0] + nu * strain[1]);
                stress[idx * 3 + 1] = factor * ((1 - nu) * strain[1] + nu * strain[0]);
                stress[idx * 3 + 2] = factor * (1 - 2 * nu) / 2 * strain[2];
            }
            ++idx;
        }
    }
}


// ����4����˹���ʵ������
void CQ8::ElementGauss(double* gp_coords)
{
    // 2x2��˹��

    int idx = 0;
    for (int i = 0; i < 2; ++i) {
        double xi = GP_Reduced[i];
        for (int j = 0; j < 2; ++j) {
            double eta = GP_Reduced[j];
            double N[8];
            ShapeFunction(xi, eta, N, nullptr, nullptr);
            double X = 0.0, Y = 0.0;
            for (int k = 0; k < 8; ++k) {
                X += N[k] * nodes_[k]->XYZ[0];
                Y += N[k] * nodes_[k]->XYZ[1];
            }
            gp_coords[idx * 2 + 0] = X;
            gp_coords[idx * 2 + 1] = Y;
            ++idx;
        }
    }
}


