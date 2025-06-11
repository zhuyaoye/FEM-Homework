#include "C3D8Material.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

// Read material data from stream Input
bool C3D8Material::Read(ifstream& Input)
{
    Input >> nset >> E >> nu;
    return true;
}

// Write material data to Stream
void C3D8Material::Write(COutputter& output)
{
    output << setw(16) << E
           << setw(16) << nu << endl;
}