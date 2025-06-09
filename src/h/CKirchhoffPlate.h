#pragma once
#include "Element.h"
#include "PlateMaterial.h"
#include <Eigen/Dense>
#include <fstream>

// Kirchhoff plate element class
class CKirchhoffPlate : public CElement {
public:
    CKirchhoffPlate();
    virtual ~CKirchhoffPlate();

    // Override base class methods
    virtual bool Read(std::ifstream& input, CMaterial* materialSets, CNode* nodeList) override;
    virtual void Write(COutputter& output) override;
    virtual void ElementStiffness(double* Matrix) override;
    virtual void ElementStress(double* stress, double* Displacement) override; // 必须定义，即使为空
    virtual void ElementLoad(double* Load, double q); // 自定义扩展函数

private:
    Eigen::Matrix<double, 3, 12> ComputeBMatrix(double xi, double eta, double a, double b);
};