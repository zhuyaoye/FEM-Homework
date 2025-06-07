/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"

#include <string>
#include <sstream>
#include "B31Material.h" // CSJ

using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	NodeListBeam = nullptr; // CSJ

	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;
	delete [] NodeListBeam; // CSJ

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::GetInstance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::GetInstance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

	// 从输入流读取一行文本
	bool hasDOF_INDEX = false;
	std::string line;
	if (std::getline(Input, line)) {
		// 创建字符串流，用于解析读取的行
		std::istringstream iss(line);
		
		// 尝试读取前四个必选参数
		if (iss >> NUMNP >> NUMEG >> NLCASE >> MODEX) {
			// 尝试读取可选的第五个参数
			if (iss >> DOF_INDEX) {
				// 成功读取了全部五个参数
				hasDOF_INDEX = true;
			} else {
				// 只有四个参数，DOF_INDEX使用默认值
				hasDOF_INDEX = false;
				DOF_INDEX = 0; // 设置默认值
				iss.clear(); // 清除错误标志
			}
			
			// 检查行中是否有多余的数据
			char extra;
			if (iss >> extra) {
				std::cerr << "警告: 行中包含多余的数据!" << std::endl;
			}
		} else {
			std::cerr << "错误: 无法解析必要的四个参数!" << std::endl;
		}
	}
	//	Read nodal point data
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//	Read load data
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

//	Read element data
	if (ReadElements())
        Output->OutputElementInfo();
    else
        return false;

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{
	if(DOF_INDEX == 0)
	{
		delete[] NodeListBeam;
	//	Read nodal point data lines
		NodeList = new CNode[NUMNP];

	//	Loop over for all nodal points
		for (unsigned int np = 0; np < NUMNP; np++)
		{
			if (!NodeList[np].Read(Input))
				return false;
		
			if (NodeList[np].NodeNumber != np + 1)
			{
				cerr << "*** Error *** Nodes must be inputted in order !" << endl
				<< "   Expected node number : " << np + 1 << endl
				<< "   Provided node number : " << NodeList[np].NodeNumber << endl;
			
				return false;
			}
		}

		return true;
	}
	else
	{
		delete[] NodeList;
	//	Read nodal point data lines
		NodeListBeam = new CBeamNode[NUMNP];

	//	Loop over for all nodal points
		for (unsigned int np = 0; np < NUMNP; np++)
		{
			if (!NodeListBeam[np].Read(Input))
				return false;
		
			if (NodeListBeam[np].NodeNumber != np + 1)
			{
				cerr << "*** Error *** Nodes must be inputted in order !" << endl
				<< "   Expected node number : " << np + 1 << endl
				<< "   Provided node number : " << NodeListBeam[np].NodeNumber << endl;
			
				return false;
			}
		}

		return true;
	}
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber() // CSJmodified
{
	NEQ = 0;
	if (DOF_INDEX == 0)
	{
		for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
		{
			for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
			{
				if (NodeList[np].bcode[dof]) 
					NodeList[np].bcode[dof] = 0;
				else
				{
					NEQ++;
					NodeList[np].bcode[dof] = NEQ;
				}
			}
		}
	}
	else
	{
		for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
		{
			for (unsigned int dof = 0; dof < CBeamNode::NDF; dof++)	// Loop over for DOFs of node np
			{
				if (NodeListBeam[np].bcode[dof]) 
					NodeListBeam[np].bcode[dof] = 0;
				else
				{
					NEQ++;
					NodeListBeam[np].bcode[dof] = NEQ;
				}
			}
		}
	}
	
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
    {
        unsigned int LL;
        Input >> LL;
        
        if (LL != lcase + 1)
        {
            cerr << "*** Error *** Load case must be inputted in order !" << endl
            << "   Expected load case : " << lcase + 1 << endl
            << "   Provided load case : " << LL << endl;
            
            return false;
        }

        LoadCases[lcase].Read(Input);
    }

	return true;
}

// Read element data
bool CDomain::ReadElements() // CSJ modified
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		if (!EleGrpList[EleGrp].Read(Input))
		{
			return false;
		}
	}
			
    return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
#ifdef _DEBUG_
    COutputter* Output = COutputter::GetInstance();
    *Output << setw(9) << "Ele = " << setw(22) << "Location Matrix" << endl;
#endif

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // Generate location matrix CSJ modified
			if(DOF_INDEX == 0)
			{
				Element.GenerateLocationMatrix();
			}
            else
			{
				Element.GenerateLocationMatrixBeam();
			}
        
		
	
#ifdef _DEBUG_
            unsigned int* LocationMatrix = Element.GetLocationMatrix();
            
            *Output << setw(9) << Ele+1;
            for (int i=0; i<Element.GetND(); i++)
                *Output << setw(5) << LocationMatrix[i];
            *Output << endl;
#endif

            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }
    
    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
    *Output << endl;
	Output->PrintColumnHeights();
#endif

}

//    Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//    and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
    //    Allocate for global force/displacement vector
    Force = new double[NEQ];
    
    //  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
    
    //    Calculate column heights
    CalculateColumnHeights();
    
    //    Calculate address of diagonal elements in banded matrix
    StiffnessMatrix->CalculateDiagnoalAddress();
    
    //    Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();
    
    COutputter* Output = COutputter::GetInstance();
    Output->OutputTotalSystemData();
}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            Element.ElementStiffness(Matrix);
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter* Output = COutputter::GetInstance();
	Output->PrintStiffnessMatrix();
#endif

}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

    clear(Force, NEQ);

//	Loop over for all concentrated loads in load case LoadCase
// CSJ modified
	if(DOF_INDEX == 0)
	{
		for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
		{
			unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
			
			if(dof) // The DOF is activated
				Force[dof - 1] += LoadData->load[lnum];
		}
	}
	else
	{
		for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
		{
			unsigned int dof = NodeListBeam[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
			
			if(dof) // The DOF is activated
				Force[dof - 1] += LoadData->load[lnum];
		}
	}
	

	return true;
}

