#include "C3D8RMaterial.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

// Read material data from stream Input
bool CC3D8RMaterial::Read(ifstream& Input)
{
    Input >> nset >> E >> nu;
    return true;
}

// Write material data to Stream
void CC3D8RMaterial::Write(COutputter& output)
{
    output << setw(16) << E
           << setw(16) << nu << endl;
}