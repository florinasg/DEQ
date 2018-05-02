/*
 * diffusion_BOUNDED_CRANK_NICHOLSON.cpp
 *
 *  Created on: 01.05.2018
 *      Author: Florian Anderl
 */

#include <vector>
#include <math.h>
#include "header.h"
#include "Defines.h"



int diffusion_ABSORBING_CRANK_NICOLSON(double a, double b)
{
	double grid_const = double(GRID_CONST);

	/*determined grid interval*/
	double interval = 0.0;
	interval = double(1/(grid_const-1));


	/*initializes grid*/
	std::vector<DIF> DENS;
	for(double i = a; i <=b ; i = i + interval )
	{
		DENS.push_back(DIF());
		DENS.back().coord = i;
	}

	return 0;
}
