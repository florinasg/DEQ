/*
 * diffusion_REFLECTIVE_CRANK_NICOLSON.cpp
 *
 *  Created on: 02.05.2018
 *      Author: Flo
 */




/*
 * diffusion_BOUNDED_CRANK_NICHOLSON.cpp
 *
 *  Created on: 01.05.2018
 *      Author: Florian Anderl
 */

#include <vector>
#include <math.h>
#include "../header.h"
#include "../Defines.h"



int diffusion_REFLECTIVE_CRANK_NICOLSON(double a, double b)
{
	double grid_const = double(GRID_CONST);

	/*determined grid interval*/
	double interval = 0.0;
	interval = double(1/(grid_const-1));


	std::ofstream file;
	file.open("CRANK_NICOLSON_IMPLICIT_"+std::to_string(GRID_CONST)+".csv");


	/*initializes grid*/
	std::vector<DIF> DENS;
	for(double i = a; i <=b ; i = i + interval )
	{
		DENS.push_back(DIF());
		DENS.back().coord = i;
	}


	DENS = dirac_delta(DENS, double(xZERO), double(SPILL_MASS), double(ALPHA));

	int T = 0;
	for(double t = 1*DELTA_T; t<= double(OBS_TIME); t = t+double(DELTA_T))
	{
		/*CALL IN MODE '2'*/
		DENS = TRI_DIAGONAL_SOLVER_REFLECTIVE(DENS,2,T);

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
