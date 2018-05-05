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
#include <vector>
#include <iostream>



template <typename T>
void print(T t)
{

  std::cout << t << " " ;
}

template<typename T, typename... Args>
void print(T t, Args... args)
{

  std::cout << t << " ";
  print(args...) ;
}



int diffusion_ABSORBING_EULER_EXPLICIT(double a, double b);

int diffusion_ABSORBING_EULER_IMPLICIT(double a, double b);

int diffusion_ABSORBING_CRANK_NICOLSON(double a, double b);

int diffusion_REFLECTIVE_EULER_EXPLICIT(double a, double b);

int REFLECTIVE_ANALYTIC(double a, double b);


typedef struct DIF{
DIF() : coord(0){};
double coord = 0;
/*history of grid point (past values + corresp. time step)*/
std::vector<std::vector<double>> conc_hist;
}DIF;



std::vector<DIF>  TRI_DIAGONAL_SOLVER_ABSORBING(std::vector<DIF> b, int mode, double time_step);
std::vector<DIF>  TRI_DIAGONAL_SOLVER_REFLECTIVE(std::vector<DIF> b_full, int mode, double time_step);
std::vector<DIF> dirac_delta(std::vector<DIF> grid, double x0, double m, double alpha);



#endif /* HEADER_H_ */
