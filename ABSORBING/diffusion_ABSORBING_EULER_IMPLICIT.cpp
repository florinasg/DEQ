/*
 * diffusion_REFLECTIVE.cpp
 *
 *  Created on: 26.04.2018
 *      Author: Florian Anderl
 */


#include <vector>
#include <math.h>
#include "../header.h"
#include "../Defines.h"



int diffusion_ABSORBING_EULER_IMPLICIT(double a, double b)
{

	std::cout << "diffusion_ABSORBING_EULER_IMPLICIT started ...\n";

	double grid_const = double(GRID_CONST);

	/*determined grid interval*/
	double interval = 0.0;
	interval = double(1/(grid_const-1));
	std::cout << "Delta X: " << interval << std::endl;

	std::ofstream file;
	file.open("EULER_IMPLICIT_"+std::to_string(GRID_CONST)+".csv");


	/*initializes grid*/
	std::vector<DIF> DENS;
	for(double i = a; i <=b ; i = i + interval )
	{
		DENS.push_back(DIF());
		DENS.back().coord = i;
		//print(i);
	}


	DENS = dirac_delta(DENS, double(xZERO), double(SPILL_MASS), double(ALPHA));


	int T = 0;
	for(double t = 1*DELTA_T; t<= double(OBS_TIME); t = t+double(DELTA_T))
	{

		/*Copies initial value at boundaries for every time step -> a*/
		DENS.at(0).conc_hist.push_back(std::vector<double>());
		DENS.at(0).conc_hist.back().push_back(T);
		DENS.at(0).conc_hist.back().push_back(double(DENS.at(0).conc_hist.at(0).at(1)));

		/*Copies initial value at boundaries for every time step -> b*/
		DENS.at(GRID_CONST-1).conc_hist.push_back(std::vector<double>());
		DENS.at(GRID_CONST-1).conc_hist.back().push_back(T);
		DENS.at(GRID_CONST-1).conc_hist.back().push_back(double(DENS.at(GRID_CONST-1).conc_hist.at(0).at(1)));


		DENS = TRI_DIAGONAL_SOLVER_ABSORBING(DENS,1,T);

		for(int i = 0; i < GRID_CONST; i++)
		{
			file << DENS.at(i).conc_hist.at(T)[1] << ",";
			file.flush();
		}
		file << "\n";
		T++;
	}

	file.close();


	return 0;
}


