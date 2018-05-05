/*
 * analytic_solutions_ABSORBING.cpp
 *
 *  Created on: 04.05.2018
 *      Author: Flo
 */



#include "../Defines.h"
#include "../header.h"


#define N 50
int REFLECTIVE_ANALYTIC(double a, double b)
{
	double interval = 1/(double(GRID_CONST)-1);

	std::ofstream file;
	file.open("REFLECTIVE_ANALYTIC_SOLUTIONS"+std::to_string(GRID_CONST)+".csv");

	std::vector<double> vals;
	std::vector<std::vector<double>> sol;

	double sum = 0;
	double eigenmode = 0;
	double eigenmode_X0 = 0;

	for(double i = a; i <= b; i = i+ interval)
	{
		vals.push_back(i);
	}

	for(double t = 1*double(DELTA_T); t<=double(OBS_TIME); t = t+double(DELTA_T))
	{

		sol.push_back(std::vector<double> ());


		for(int i = a; i <= b;i = i+ interval)
		{


			for (int j = 0; j <= N; j++)
			{
				if(j == 0)
				{
					eigenmode = 0;
					eigenmode_X0 = 0;
				}
				else
				{
					eigenmode = sqrt(2/b) * sin(j*M_PI*(i/b));
					eigenmode_X0 = sqrt(2/b) * sin(j*M_PI*(xZERO*interval/b));
				}
				sum = sum + exp(-pow(((M_PI*j)/b),2)*double(DIFFUSION_CONST)*double(DELTA_T))*eigenmode*eigenmode_X0;
			}

			sol.back().push_back(double());
			sol.back().back() =

					double(SPILL_MASS) * sum;



			sum = 0;

			file << sol.back().back() << "," ;

		}

		file << "\n";

	}


	file.close();
	return 0;

}

