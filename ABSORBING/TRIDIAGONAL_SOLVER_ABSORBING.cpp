/*
 * TRIDIAGONAL_SOLVER.cpp
 *
 *  Created on: 02.05.2018
 *      Author: Flo
 */


#include "../header.h"
#include "../Defines.h"
#include <armadillo>
#include <vector>

/*IMPLEMENTS TRIDIAGONAL SOLVER ->
 * mode: 0 -> explicit
 * 1 -> implicit*/
std::vector<DIF>  TRI_DIAGONAL_SOLVER_ABSORBING(std::vector<DIF> b_full, int mode, double time_step)
{

	/*excludes first and last element -> boundary condition*/
	std::vector<DIF> b(b_full.begin()+1, b_full.begin()+GRID_CONST-1);
	//print("b:size()",b.size(),"b.back().coord", b.back().coord);
	std::vector<DIF> b_return(GRID_CONST);



	double interval = b[1].coord - b[0].coord;



	/*EULER IMPLICIT*/
	if(mode == 1)
	{
		/*Dimensions implement boundary conditions*/
		arma::mat A(GRID_CONST-2, GRID_CONST-2);
		arma::vec B(GRID_CONST-2);
		//print("B.size()",B.size());

		arma::vec X;


		A.zeros();


		/*constructs vector B*/
		for(int i = 0; i < B.size(); i++)
		{
			B(i) = b.at(i).conc_hist.back().at(1);
			//print(i, b.at(i).conc_hist.back().at(1));

		}


		/*Construct A MATRIX*/


		double alpha = double(DIFFUSION_CONST) * (double(DELTA_T)/double(pow(interval,2)));

		for(int i = 0; i < double(GRID_CONST)-2; i++)
		{
			if(i == 0)
			{
				A(i,i) = 1+2*alpha;
				A(i,i+1) = - alpha;
			}

			else if(i == double(GRID_CONST)-3)
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

		//A.save("A_EULER_IMPLICIT.csv",arma::csv_ascii);
		//B.save("B_EULER_IMPLICIT.csv", arma::csv_ascii);
		X = arma::solve(A, B,arma::solve_opts::fast);



		for(int i = 0; i < b.size(); i++)
		{
			b.at(i).conc_hist.push_back(std::vector<double>());
			b.at(i).conc_hist.back().push_back(time_step+1);
			b.at(i).conc_hist.back().push_back(double(X(i)));
			print(i);

		}
	}

	/*CRANK-NICOLSON*/
	else if(mode == 2)
	{

		arma::mat A(GRID_CONST-2, GRID_CONST-2);
		arma::mat B1(GRID_CONST-2,GRID_CONST-2);
		arma::vec B2(GRID_CONST-2);

		arma::mat B(GRID_CONST-2,GRID_CONST-2);

		arma::vec X;
		A.zeros();


		double alpha = double(DIFFUSION_CONST) * (double(DELTA_T)/double(pow(interval,2)));

		/*CONSTRUCT A MATRIX*/
		for(int i = 0; i < double(GRID_CONST)-2; i++)
		{
			if(i == 0)
			{
				A(i,i) = 1+alpha;
				A(i,i+1) = - (alpha/2);
			}

			else if(i == double(GRID_CONST)-3)
			{
				A(i,i) = 1+alpha;
				A(i,i-1) = - (alpha/2);
			}

			else
			{
				A(i,i) = 1+alpha;
				A(i,i+1) = - (alpha/2);
				A(i,i-1) = - (alpha/2);
			}


		}



		/*CONSTRUCT B1-MATRIX*/
		for(int i = 0; i < double(GRID_CONST)-2; i++)
		{
			if(i == 0)
			{
				B1(i,i) = 1-alpha;
				B1(i,i+1) = (alpha/2);
			}

			else if(i == double(GRID_CONST)-3)
			{
				B1(i,i) = 1-alpha;
				B1(i,i-1) = (alpha/2);
			}

			else
			{
				B1(i,i) = 1-alpha;
				B1(i,i+1) = (alpha/2);
				B1(i,i-1) =  (alpha/2);
			}


		}


		/*constructs vector B2*/
		for(int i = 0; i < B2.size(); i++)
		{
			B2(i) = b.at(i).conc_hist.back().at(1);

		}

		/*TRANSFORMS B-MATRICES INTO VECTOR*/
		B = B1 * B2;


		X = arma::solve(A, B, arma::solve_opts::fast);


		print("\nb: ",b.size());
		for(int i = 0; i < b.size(); i++)
		{
			b.at(i).conc_hist.push_back(std::vector<double>());
			b.at(i).conc_hist.back().push_back(time_step+1);
			b.at(i).conc_hist.back().push_back(double(X(i)));

		}
	}

	print("\nb_full: ", b_full.size());
	print("\nb_return: ", b_return.size());

	b_return[0] = b_full[0];

	for(int i = 1; i < b_return.size()-1; i++)
	{
		b_return[i] = b[i-1];
	}

	b_return[b_return.size()-1] = b_full.back();


	return b_return;
}


