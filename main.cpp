/*
 * main.cpp
 *
 *  Created on: 21.04.2018
 *      Author: Florian Anderl
 */

#include "header.h"


int main(int argc, char * argv[])
{
	/*args = [a,b, diffusion_constant, x0(center of dirac-sist) factor*grid_points(e.g. 0.5 * 1000), mass spilled]*/
	diffusion_BOUNDED_EULER_EXPLICIT(0,1);

	return 0;
}
