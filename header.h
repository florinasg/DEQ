/*
 * header.h
 *
 *  Created on: 26.04.2018
 *      Author: Florian Anderl
 */

#ifndef HEADER_H_
#define HEADER_H_

#include <armadillo>
#include <cmath>
#include <ctgmath>
#include <math.h>
#include <fstream>


int diffusion_BOUNDED_EULER_EXPLICIT(double a, double b);

int diffusion_BOUNDED_EULER_IMPLICIT(double a, double b);

int diffusion_BOUNDED_EULER_CRANK_NICOLSON(double a, double b);



typedef struct DIF{
DIF() : coord(0){};
double coord = 0;
/*history of grid point (past values + corresp. time step)*/
std::vector<std::vector<double>> conc_hist;
}DIF;




std::vector<DIF> diract_delta(std::vector<DIF> grid, double x0, double m, double alpha);



#endif /* HEADER_H_ */
