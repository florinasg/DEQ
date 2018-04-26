/*
 * diffusion_REFLECTIVE.cpp
 *
 *  Created on: 26.04.2018
 *      Author: Florian Anderl
 */

#include <vector>
#include <math.h>
#include "header.h"



/*Defines number of grid points*/
#define GRID_CONST 1000
#define x0 0.5
#define ALPHA 0.0000000000000001
#define SPILL_MASS 10
#define DIFFUSION_CONST 0.5
#define OBS_TIME 10000


/*EXPLICIT -> NO LINEAR SYSTEM*/

int diffusion_BOUNDED_EULER_EXPLICIT(double a, double b)
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
	/*UNIT TEST*/
	//	for(int i = 0; i < DENS.size(); i++)
	//	{
	//	std::cout << DENS.at(i).coord << std::endl;
	//	}


	/*implements application of initial conditon */
	DENS = diract_delta(DENS, double(x0)*grid_const, double(SPILL_MASS), double(ALPHA));
	/*UNIT TEST*/
	//	for(int i = 0; i < DENS.size(); i++)
	//	{
	//	std::cout << DENS.at(i).density << std::endl;
	//	}


	/*Implementation & Iteration of EXPLICIT EULER SCHEME*/
	for(double t = 1; t<= double(OBS_TIME); t = t+1)
	{
		/*Boundrary Condition*/
		for(int i = 1; i < DENS.size()-1; i++)
		{
			if(i == 1)
			{
				/*adds element in density history of DENS(i)*/
				DENS.at(i).conc_hist.push_back(std::vector<double>());
				/*writes time step into element*/
				DENS.at(i).conc_hist.back().push_back(t);
				/*writes updated value of density in element*/
				DENS.at(i).conc_hist.back().push_back(double());
				DENS.at(i).conc_hist.back()[2] =

						/*EULER EXPLICIT*/
						DENS.at(i).conc_hist.at(t-1)[2] +
						double(DIFFUSION_CONST)* (1/pow(interval,2)) * (DENS.at(i+1).conc_hist.at(t-1)[2]
																									   - 2*DENS.at(i).conc_hist.at(t-1)[2]);



			}
			else if(i == DENS.size()-2)
			{

				/*adds element in density history of DENS(i)*/
				DENS.at(i).conc_hist.push_back(std::vector<double>());
				/*writes time step into element*/
				DENS.at(i).conc_hist.back().push_back(t);
				/*writes updated value of density in element*/
				DENS.at(i).conc_hist.back().push_back(double());
				DENS.at(i).conc_hist.back()[2] =

						/*EULER EXPLICIT*/
						DENS.at(i).conc_hist.at(t-1)[2] +
						double(DIFFUSION_CONST)* (1/pow(interval,2)) * (- 2*DENS.at(i).conc_hist.at(t-1)[2] + DENS.at(i-1).conc_hist.at(t-1)[2]);

			}


			else
			{
				/*adds element in density history of DENS(i)*/
				DENS.at(i).conc_hist.push_back(std::vector<double>());
				/*writes time step into element*/
				DENS.at(i).conc_hist.back().push_back(t);
				/*writes updated value of density in element*/
				DENS.at(i).conc_hist.back().push_back(double());
				DENS.at(i).conc_hist.back()[2] =

						/*EULER EXPLICIT*/
						DENS.at(i).conc_hist.at(t-1)[2] +
						double(DIFFUSION_CONST)* (1/pow(interval,2)) * (DENS.at(i+1).conc_hist.at(t-1)[2]
																									   - 2*DENS.at(i).conc_hist.at(t-1)[2]
																																		+ DENS.at(i-1).conc_hist.at(t-1)[2]);

			}

		}
	}

	/*writes DESNITIES INTO FILE*/

	std::ofstream file;
	file.open("EULER_EXPLICIT.csv");

	for(int i = 0; i < DENS.size(); i++)
	{
		for(int t  = 0; t< double(OBS_TIME); t++)
		{
			file << DENS.at(i).conc_hist.at(t)[2] << "\n";
		}

		file << ",";
	}


	file.close();

			return 0;
}








int diffusion_BOUNDED_EULER_IMPLICIT(double a, double b)
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




int diffusion_BOUNDED_EULER_CRANK_NICOLSON(double a, double b)
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

