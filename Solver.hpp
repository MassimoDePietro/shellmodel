/*
* Copyright 2017 Massimo De Pietro
*
* This file is free software: you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the Free Software
* Foundation, either version 3 of the License, or (at your option) any later
* version.
*
* This file is distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
* A PARTICULAR PURPOSE. See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with
* this file. If not, see http://www.gnu.org/licenses/.
*
*
* DESCRIPTION:
* This scientific code solves the equations of motion for the "SABRA"
* shell model (see https://doi.org/10.1103/PhysRevE.58.1811).
*
*
* DEPENDENCIES:
* None. Just a compiler supporting c++11.
*
* VERSION:
* v3
*
*/


/// \file Solver.hpp
/// \brief This file contains the declaration of the Solver class, the Solution 
/// class and structures for the parameters. The Solution class stores just the 
/// values of time and a velocity field. The solver class is the object that 
/// carry on the simulation and the calculation of the values of the observables.



#pragma once


#include <complex>
#include <vector>

#include "Exception_with_message.hpp"
#include "File.hpp"
#include "Observables.hpp"




struct Physical_parameters {
	unsigned int N;
	std::complex<double> F0, F1;
	double nu;
	// Hardcoded shell model parameters
	const double lam = 2.0;	//Lambda
	const double K0 = 1.0;
	const double A = +1.0;
	const double B = -0.5;
	const double C = +0.5;
};

struct Control_parameters {
	double dt;
	unsigned long long ntimesteps;
	unsigned long long calcsteps;
	unsigned long long dumpsteps;
};




class Solution 
{
private:
	const unsigned int default_initialization_n_max = 6;
	const double initialization_K0 = 1.0;
	const double initialization_lambda = 2.0;

public:
	Solution(unsigned int N);											///< Init all u's to 0
	Solution(unsigned int N, unsigned long int random_seed);			///< Init all u's to random
	Solution(unsigned int N, std::string filename);						///< Init reading the u's from file
	~Solution();

	double t;
	std::vector< std::complex<double> > u;
};


class Solver
{
public:
	// Constructor/destructor
	Solver() {};
	~Solver();
	
	// Physical attributes
	unsigned int N;
	std::complex<double> F0, F1;
	double nu;
	const double lam = 2.0;	//Lambda
	const double K0 = 1.0;
	const double A = +1.0;
	const double B = -0.5;
	const double C = +0.5;
	std::vector<double> k;

	// Solution at time t
	Solution* s;			// Contains the solution at time t

	// Observables
	observables current;    // Observables, instantaneous values
	observables stat;       // Observables, time averages
	
	// Internal copy of the observables tables
	std::vector<f_of_t_table_record> table_functions_of_t;
	std::vector<f_of_t_n_table_record> table_functions_of_t_n;

    // Attributes for integration control 
    double dt;
    long long ntimesteps, calcsteps, dumpsteps;
    long long statcnt;
    long long stepcnt;

	// Temporary memory for computations
	std::vector<double> k2;
	std::vector<std::complex<double>> du1, du2, du3, du4;     // For intermediate runge kutta steps
	std::vector<std::complex<double>> u2, u3, u4;			  // For intermediate runge kutta steps
    
	// Methods
	void init(Physical_parameters params, Control_parameters control, const std::vector<f_of_t_table_record>& l_table_functions_of_t,
		const std::vector<f_of_t_n_table_record>& l_table_functions_of_t_n);
	void assign_solution(Solution s0);
	void assign_statistics(observables& stat0);

	/// nlt() with no arguments calculates the non linear term using the current solution
	void nlt();
	/// nlt(...) with the arguments calculates the non linear term using the solution passed as argument. Consider it as a static method.
	void nlt(const std::vector<std::complex<double>>& u, const std::vector<double>& k, std::vector<std::complex<double>>& du);

	/// Right hand side of the differential equation: nlt + dissipation + forcing
	void rhs(const std::vector<std::complex<double>>& u, const std::vector<double>& k, const std::vector<double>& k2, std::complex<double> F0, std::complex<double> F1, double nu, std::vector<std::complex<double>>& du);

    void advance_one_step_rk4();

    void calculate_observables();
};

