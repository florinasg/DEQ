/*
 * diffusion_REFLECTIVE_EULER_EXPLICIT.cpp
 *
 *  Created on: 02.05.2018
 *      Author: Flo
 */




/*
 * diffusion_BOUNDEN_EULER_EXPLICIT.cpp
 *
 *  Created on: 01.05.2018
 *      Author: Florian Anderl
 */

#include <vector>
#include <math.h>
#include "../header.h"
#include "../Defines.h"






/*EXPLICIT -> NO LINEAR SYSTEM*/

int diffusion_ABSORBING_EULER_EXPLICIT(double a, double b)
{

	std::cout << "EULER_EXPLICIT SCHEME WITH ABSORBING BOUNDRARIER... " << std::endl;

	double grid_const = double(GRID_CONST);

	/*determined grid interval*/
	double interval = 0.0;
	interval = double(1/(grid_const-1));
	std::cout << "Delta X: " << interval << std::endl;


	/*defines time-step*/
	double deltaT = double(DELTA_T);





	/*initializes grid*/
	std::vector<DIF> DENS;
	for(double i = a; i <= b ; i = i + interval )
	{
		DENS.push_back(DIF());
		DENS.back().coord = i;
		std::cout << i << std::endl;
	}



	/*implements application of initial conditon */
	DENS = dirac_delta(DENS, double(xZERO), double(SPILL_MASS), double(ALPHA));



	/*Implementation & Iteration of EXPLICIT EULER SCHEME*/


	std::ofstream file;
	file.open("EULER_EXPLICIT_"+std::to_string(GRID_CONST)+".csv");


	for(int i = 1; i < DENS.size()-1; i++)
	{
		file << DENS.at(i).conc_hist.at(0)[1] << ",";
	}

	file << "\n";

	int T = 1;
	for(double t = 1*DELTA_T; t<= double(OBS_TIME); t = t+double(DELTA_T))
	{
		/*Copies initial value at boundaries for every time step -> a*/
		DENS.at(0).conc_hist.push_back(std::vector<double>());
		DENS.at(0).conc_hist.back().push_back(T);
		DENS.at(0).conc_hist.back().push_back(double(DENS.at(0).conc_hist.at(0).at(1)));

		/*Copies initial value at boundaries for every time step -> b*/
		DENS.at(GRID_CONST-1).conc_hist.push_back(std::vector<double>());
		DENS.at(GRID_CONST-1).conc_hist.back().push_back(T);
		DENS.at(GRID_CONST-1).conc_hist.back().push_back(double(DENS.at(GRID_CONST-1).conc_hist.at(0).at(1)));


		/*Boundary Condition implemented in boundaries of for-loop*/
		for(int i = 1; i < DENS.size()-1; i++)
		{
			if(i == 1)
			{
				/*adds element in density history of DENS(i)*/
				DENS.at(i).conc_hist.push_back(std::vector<double>());
				/*writes time step into element*/
				DENS.at(i).conc_hist.back().push_back(T);
				/*writes updated value of density in element*/
				DENS.at(i).conc_hist.back().push_back(double());
				double nominator = pow(interval,2);
				DENS.at(i).conc_hist.back()[1] =

						/*EULER EXPLICIT*/
						DENS.at(i).conc_hist.at(T-1)[1] +
						double(DIFFUSION_CONST)* (deltaT/nominator) * (DENS.at(i+1).conc_hist.at(T-1)[1]
																									  - 2*DENS.at(i).conc_hist.at(T-1)[1]);

			}

			else if(i == DENS.size()-2)
			{

				/*adds element in density history of DENS(i)*/
				DENS.at(i).conc_hist.push_back(std::vector<double>());
				/*writes time step into element*/
				DENS.at(i).conc_hist.back().push_back(T);
				/*writes updated value of density in element*/
				DENS.at(i).conc_hist.back().push_back(double());
				double nominator = pow(interval,2);
				DENS.at(i).conc_hist.back()[1] =

						/*EULER EXPLICIT*/
						DENS.at(i).conc_hist.at(T-1)[1] +
						double(DIFFUSION_CONST)* (deltaT/nominator) * (- 2*DENS.at(i).conc_hist.at(T-1)[1] + DENS.at(i-1).conc_hist.at(T-1)[1]);

			}


			else
			{
				/*adds element in density history of DENS(i)*/
				DENS.at(i).conc_hist.push_back(std::vector<double>());
				/*writes time step into element*/
				DENS.at(i).conc_hist.back().push_back(T);
				/*writes updated value of density in element*/
				DENS.at(i).conc_hist.back().push_back(double());
				double nominator = pow(interval,2);
				DENS.at(i).conc_hist.back()[1] =

						/*EULER EXPLICIT*/
						DENS.at(i).conc_hist.at(T-1)[1] +
						double(DIFFUSION_CONST)* (deltaT/nominator) *
						(DENS.at(i+1).conc_hist.at(T-1)[1]

														- 2*DENS.at(i).conc_hist.at(T-1)[1]

																						 + DENS.at(i-1).conc_hist.at(T-1)[1]);

			}


			/*save on Harddrive*/
			file << DENS.at(i).conc_hist.at(T)[1] << ",";
			file.flush();




		}

		file << "\n";

		T++;
	}

	/*writes DESNITIES INTO FILE*/

	//	std::ofstream file;
	//	file.open("EULER_EXPLICIT_"+std::to_string(GRID_CONST)+".csv");
	//	std::string prefix = "";
	//	T = 0;
	//
	//	for(double t  = 0; t< double(OBS_TIME); t = t + double(DELTA_T))
	//	{
	//
	//		for(int i = 0; i < DENS.size(); i++)
	//		{
	//			file << DENS.at(i).conc_hist.at(T)[1] << ",";
	//		}
	//
	//		file << "\n";
	//		T++;
	//	}
	//	file.close();


	file.close();
	return 0;
}



