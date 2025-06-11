/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "ElementGroup.h"
#include "Domain.h"

CNode* CElementGroup::NodeList_ = nullptr;

//! Constructor
CElementGroup::CElementGroup()
{
    if (!NodeList_)
    {
        CDomain* FEMData = CDomain::GetInstance();
        NodeList_ = FEMData->GetNodeList();
    }
    
    ElementType_ = ElementTypes::UNDEFINED;
    
    NUME_ = 0;
    ElementList_ = nullptr;
    
    NUMMAT_ = 0;
    MaterialList_ = nullptr;
}

//! Deconstructor
CElementGroup::~CElementGroup()
{
    if (ElementList_)
        delete [] ElementList_;
    
    if (MaterialList_)
        delete [] MaterialList_;
}

//! operator []
//! For the sake of efficiency, the index bounds are not checked
CElement& CElementGroup::operator[](unsigned int i)
{
    return *(CElement*)((std::size_t)(ElementList_) + i*ElementSize_);
}

//! Return index-th material in this element group
CMaterial& CElementGroup::GetMaterial(unsigned int i)
{
    return *(CMaterial*)((std::size_t)(MaterialList_) + i*MaterialSize_);
}

//! Calculate the size of the derived element and material class
void CElementGroup::CalculateMemberSize()
{
    switch (ElementType_)
    {
        case ElementTypes::UNDEFINED:
            std::cerr << "Setting element type to UNDEFINED." << std::endl;
            exit(5);
        case ElementTypes::Bar:
            ElementSize_ = sizeof(CBar);
            MaterialSize_ = sizeof(CBarMaterial);
            break;
        case ElementTypes::Plate:
            ElementSize_ = sizeof(CKirchhoffPlate);
            MaterialSize_ = sizeof(CPlateMaterial);
            break;
        case ElementTypes::Beam: //CSJ
            ElementSize_ = sizeof(CB31); //CSJ
            MaterialSize_ = sizeof(CB31Material); //CSJ
            break; //CSJ
        case ElementTypes::H8:
            ElementSize_ = sizeof(C3D8);
            MaterialSize_ = sizeof(C3D8Material);
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::CalculateMemberSize." << std::endl;
            exit(5);
            break;
    }
}

//! Allocate array of derived elements
void CElementGroup::AllocateElements(std::size_t size)
{
    switch(ElementType_)
    {
        case ElementTypes::Bar:
            ElementList_ = new CBar[size];
            break;
        case ElementTypes::Plate:
            ElementList_ = new CKirchhoffPlate[size];
            break;
        case ElementTypes::Beam: //CSJ
            ElementList_ = new CB31[size]; //CSJ
            break; //CSJ
        case ElementTypes::H8:
            ElementList_ = new C3D8[size];
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::AllocateElement." << std::endl;
            exit(5);
    }
}

//! Allocate array of derived materials
void CElementGroup::AllocateMaterials(std::size_t size)
{
    switch(ElementType_)
    {
        case ElementTypes::Bar:
            MaterialList_ = new CBarMaterial[size];
            break;
        case ElementTypes::Plate:
            MaterialList_ = new CPlateMaterial[size];
            break;
        case ElementTypes::Beam: //CSJ
            MaterialList_ = new CB31Material[size]; //CSJ
            break; //CSJ
        case ElementTypes::H8:
            MaterialList_ = new C3D8Material[size];
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::AllocateMaterial." << std::endl;
            exit(5);
    }
}

//! Read element group data from stream Input
bool CElementGroup::Read(ifstream& Input)
{
    Input >> (int&)ElementType_ >> NUME_ >> NUMMAT_;
    
    CalculateMemberSize();

//  Read material/section property lines
    AllocateMaterials(NUMMAT_);
    
//  Loop over for all material property sets in this element group
    for (unsigned int mset = 0; mset < NUMMAT_; mset++)
    {
        GetMaterial(mset).Read(Input);
  
        if (GetMaterial(mset).nset != mset + 1)
        {
            cerr << "*** Error *** Material sets must be inputted in order !" << endl
            << "    Expected set : " << mset + 1 << endl
            << "    Provided set : " << GetMaterial(mset).nset << endl;
        
            return false;
        }
    }

//  Read element data lines
    AllocateElements(NUME_);
    
//  Loop over for all elements in this element group
    for (unsigned int Ele = 0; Ele < NUME_; Ele++)
    {
        unsigned int N;
        
        Input >> N;    // element number
        
        if (N != Ele + 1)
        {
            cerr << "*** Error *** Elements must be inputted in order !" << endl
            << "    Expected element : " << Ele + 1 << endl
            << "    Provided element : " << N << endl;
            
            return false;
        }

        if (!(*this)[Ele].Read(Input, MaterialList_, NodeList_))
            return false;
    }

    return true;
}
