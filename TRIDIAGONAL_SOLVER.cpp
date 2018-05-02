/*
 * TRIDIAGONAL_SOLVER.cpp
 *
 *  Created on: 02.05.2018
 *      Author: Flo
 */


#include "header.h"
#include "Defines.h"
#include <armadillo>

/*IMPLEMENTS TRIDIAGONAL SOLVER ->
 * mode: 0 -> explicit
 * 1 -> implicit*/
std::vector<DIF>  TRI_DIAGONAL_SOLVER(std::vector<DIF> b, int mode, double time_step)
{


	double interval = b[1].coord - b[0].coord;

	arma::mat A(GRID_CONST, GRID_CONST);
	arma::vec B(GRID_CONST,1);

	arma::vec X;


	A.zeros();

	/*EULER IMPLICIT*/
	if(mode == 1)
	{


		/*constructs vector B*/
		for(int i = 0; i < B.size(); i++)
		{
			B(i,1) = b.at(i).conc_hist.back().at(1);
		}


		/*Construct A MATRIX*/


		double alpha = double(DIFFUSION_CONST) * (double(DELTA_T)/double(pow(interval,2)));

		for(int i = 1; i <= A.size(); i++)
		{
			if(i == 1)
			{
				A(i,i) = 1+2*alpha;
				A(i,i+1) = - alpha;
			}

			else if(i == A.size())
			{
				A(i,i) = 1+2*alpha;
				A(i,i-1) = - alpha;
			}

			else
			{
				A(i,i) = 1+2*alpha;
				A(i,i+1) = - alpha;
				A(i,i-1) = - alpha;
			}
		}



	}


	X = arma::solve(A, B);



	for(int i = 0; i < b.size(); i++)
	{
		b.at(i).conc_hist.push_back(std::vector<double>());
		b.at(i).conc_hist.back().push_back(time_step+1);
		b.at(i).conc_hist.back().push_back(double(X(i,1)));

	}




	return b;
}


