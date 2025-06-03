/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"
#include <Eigen/Dense>
using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

public:

//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output) = 0;

	// 用于 Kirchhoff 板单元的函数：声明为虚函数 + 默认实现（抛异常）
    virtual Eigen::Matrix3d GetFlexuralMatrix() const {
        throw std::runtime_error("GetFlexuralMatrix() not supported by this material.");
    }

    virtual double GetThickness() const {
        throw std::runtime_error("GetThickness() not supported by this material.");
    }

    virtual Eigen::Matrix3d GetstressMatrix() const {
        throw std::runtime_error("GetstressMatrix() not supported by this material.");
    }
};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};
