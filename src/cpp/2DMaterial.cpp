

#include "2DMaterial.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool C2DMaterial::Read(ifstream& Input)
{
	Input >> nset;	

	Input >> E >> nu >> t >> type;

	return true;
}


//	Write material data to Stream
void C2DMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16)<<nu<<setw(16)<<t<<setw(16) << type << endl;
}
