/*
 * dirac_delta.cpp
 *
 *  Created on: 26.04.2018
 *      Author: Florian Anderl
 */


#include "header.h"
#include "Defines.h"

/*implementation of Dirac-delta distribution
 * x0 = center of distribution
 * returns grid 'treated with dirac'*/

/*args=[grid,distribution center, spill-mass, alpha (see dirac)]*/
std::vector<DIF> dirac_delta(std::vector<DIF> grid, double x0, double m, double alpha)
{

	double interval = grid[1].coord - grid[0].coord;


	/*TODO: OLD IMPLEMENTATION OF DIRAC*/
	/*applies dirac delta distribution */
	//	for(int i = 0; i< grid.size(); i++)
	//	{
	//		grid.at(i).conc_hist.push_back(std::vector<double>());
	//		grid.at(i).conc_hist.back().push_back(0);
	//		grid.at(i).conc_hist.back().push_back(double());
	//		grid.at(i).conc_hist.back()[1] =
	//				m * (1/(alpha*sqrt(M_PI))) * exp(-
	//						pow((grid.at(i).coord-(grid.at(x0).coord)),2) /
	//						pow(alpha,2)) ;
	//	}

	/*NEW IMPLEMENTATION GOES HERE*/
	for(int i = 0; i< grid.size(); i++)
	{
		if(i == x0)
		{
			grid.at(i).conc_hist.push_back(std::vector<double>());
			grid.at(i).conc_hist.back().push_back(double(1*(DELTA_T)));
			grid.at(i).conc_hist.back().push_back(double());
			grid.at(i).conc_hist.back()[1] = m * 1/interval;
		}

		else
		{
			grid.at(i).conc_hist.push_back(std::vector<double>());
			grid.at(i).conc_hist.back().push_back(0);
			grid.at(i).conc_hist.back().push_back(double());
			grid.at(i).conc_hist.back()[1] = m * 0;
		}

	}


	return grid;
}
