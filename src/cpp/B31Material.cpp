#include "B31Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Read material data from stream Input
bool CB31Material::Read(ifstream& Input)
{
	// Input >> nset >> E >> nu >> A;
	// G= E / (2 * (1 + nu));
	// // 尝试读取 Iy, Iz, J，如果失败，则自动计算
	// if (!(Input >> Iy >> Iz >> J)) {
	// 	Input.clear(); // 清除状态
	// 	double H, B;
	// 	if (Input >> H >> B) {
	// 		// 以矩形截面为例自动计算 Iy, Iz, J
	// 		Iy = B * pow(H, 3) / 12.0;
	// 		Iz = H * pow(B, 3) / 12.0;
	// 		//J = (B * pow(H, 3)) / 3.0; // 或者使用更合适的近似公式
	// 		double beta = B / H;
	// 		J = (B * pow(H, 3)) / 3.0 * (1.0 - 0.21 * beta * (1 - pow(beta, 4) / 12.0));
	// 	} else {
	// 		cerr << "Error: Missing Iy, Iz, J or H, B for section properties." << endl;
	// 		return false;
	// 	}
	// }
	double H, B, t_H, t_B;
	Input >> nset >> E >> nu >> A >> H >> B >> t_H >> t_B;
	G= E / (2 * (1 + nu));
	double H_in, B_in;
	H_in = H - 2*t_H;
	B_in = B - 2*t_B;
	// 以矩形截面为例自动计算 Iy, Iz, J
	Iy = B * pow(H, 3) / 12.0 - B_in * pow(H_in, 3) / 12.0;
	Iz = H * pow(B, 3) / 12.0 - H_in * pow(B_in, 3) / 12.0;
	//J = (B * pow(H, 3)) / 3.0; // 或者使用更合适的近似公式
	double beta = B / H;
	J = (B * pow(H, 3)) / 3.0 * (1.0 - 0.21 * beta * (1 - pow(beta, 4) / 12.0)) - (B_in * pow(H_in, 3)) / 3.0 * (1.0 - 0.21 * beta * (1 - pow(beta, 4) / 12.0));

	double BHArea;
	BHArea = B*H;
    if (BHArea != A)
    {
        cerr << "B and H do not fit with A." << endl;
        exit(EXIT_FAILURE);
    }
	return true;
}


//	Write material data to Stream
void CB31Material::Write(COutputter& output)
{
	output << setw(16) << E
	       << setw(16) << nu
	       << setw(16) << A
	       << setw(16) << Iy
	       << setw(16) << Iz
	       << setw(16) << J << endl;
}