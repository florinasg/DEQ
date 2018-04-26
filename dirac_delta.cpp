/*
 * dirac_delta.cpp
 *
 *  Created on: 26.04.2018
 *      Author: Florian Anderl
 */


#include "header.h"

/*implementation of Dirac-delta distribution
 * x0 = center of distribution
 * returns grid 'treated with dirac'*/

/*args=[grid,distribution center, spill-mass, alpha (see dirac)]*/
std::vector<DIF> diract_delta(std::vector<DIF> grid, double x0, double m, double alpha)
{


	/*dummy grid as array -> NOT USED CURRENTLY*/
	std::vector<double> grid_DUM;
	grid_DUM.resize(grid.size()); /*TODO: TEST*/


	/*applies dirac delta distribution */
	for(int i = 0; i< grid_DUM.size(); i++)
	{
		grid.at(i).conc_hist.push_back(std::vector<double>());
		grid.at(i).conc_hist.back().push_back(0);
		grid.at(i).conc_hist.back().push_back(double());
		grid.at(i).conc_hist.back().at(1) =
		 m * (1/(alpha*sqrt(M_PI))) * exp(-
				pow((grid.at(i).coord-(grid.at(x0).coord)),2) /
				pow(alpha,2)) ;
	}
	/*TODO: VERIFY DIRAC*/



	return grid;
}
