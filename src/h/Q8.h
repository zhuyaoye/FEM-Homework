
#pragma once

#include "Element.h"


class CQ8 : public CElement
{
public:

    //!	Constructor
    CQ8();

    //!	Desconstructor
    ~CQ8();

    //!	Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

    //!	Write element data to stream
    virtual void Write(COutputter& output);

    //!	Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

    //!	Calculate element stress
    virtual void ElementStress(double* stress, double* Displacement);
    void ElementGauss(double* gp_coords); // 计算4个高斯点的实际坐标

private:

    void ShapeFunction(double xi, double eta, double* N, double* dNdxi, double* dNdeta); // 形函数
    void Jacobian(double xi, double eta, double& detJ, double* dNdx, double* dNdy);     // 雅可比矩阵
    void GaussIntegration(double* Matrix);                                              // 高斯积分
   
};

   