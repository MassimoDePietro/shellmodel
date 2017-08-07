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


/// \file Solver.cpp
/// \brief This file contains the definition of the Solver class and the the 
/// Solution class.


#include "Solver.hpp"
#include "check.hpp"

//#include <array>
//#include <cmath>
#include <cstdio>
#include <cstdlib>
//#include <cstring>
#include <iostream>
//#include <map>
//#include <stdexcept>
//#include <vector>

using namespace std;


complex<double> C0 = complex<double>(0., 0.);
complex<double> I  = complex<double>(0., 1.);

const double PI = 3.14159265358979323846264338327;
const double TWOPI = 2.*PI;



Solver::~Solver()
{
	delete s;
}

void Solver::nlt()
{
	nlt(s->u, k, du1);
}

void Solver::nlt(const std::vector<std::complex<double>>& u, const std::vector<double>& k, std::vector<std::complex<double>>& du)
{
    // u is the vector with the velocities
    // du is the vector where non linear term will be written
    
    complex<double> aa,bb,cc;
	const double coeffa = lam * A;
	const double coeffb = B;
	const double coeffc = C / lam;
    
    const int N = u.size();
    for (int j=0; j<N; j++)
    {
        if ( (j+2) < N ) aa = conj(u[j+1])*u[j+2];
        else aa = C0;   // C0 = complex<double>(0,0)
        
        if (( (j+1) < N ) && ( (j-1) >= 0 ) ) bb = conj(u[j-1])*u[j+1];
        else bb = C0;
        
        if ( (j-2) >= 0 ) cc =  u[j-2]*u[j-1];
        else cc = C0;
        
        du[j]=I*k[j]*(coeffa*aa + coeffb*bb + coeffc*cc);
    }
}


void Solver::rhs(const vector<complex<double>>& u, const vector<double>& k, const vector<double>& k2, complex<double> F0, complex<double> F1, double nu, vector<complex<double>>& du)
{
	// u is the vector with the velocities
	// du is the vector where the derivative will be written

	// NLT
	nlt(u, k, du);
	
	// Forcing
	du[0] += F0;
	du[1] += F1;

	// Viscosity
	for (int j = 0; j < u.size(); j++) du[j] -= nu*k2[j]*u[j];
}



void Solver::advance_one_step_rk4()
{
	// Runge kutta 4
	const double dt2 = dt / 2.0;
	const double dt3 = dt / 3.0;
	const double dt6 = dt / 6.0;

	stepcnt++;

	// K1
	rhs(s->u, k, k2, F0, F1, nu, du1);

	// K2
	for (int j = 0; j < N; j++) {
		u2[j] = s->u[j] + dt2*du1[j];
	}
	rhs(u2, k, k2, F0, F1, nu, du2);

	// K3
	for (int j = 0; j < N; j++) {
		u3[j] = s->u[j] + dt2*du2[j];
	}
	rhs(u3, k, k2, F0, F1, nu, du3);

	// K4
	for (int j = 0; j < N; j++) {
		u4[j] = s->u[j] + dt*du3[j];
	}
	rhs(u4, k, k2, F0, F1, nu, du4);

	// Advance solution
	for (int j = 0; j < N; j++) {
		s->u[j] += dt6*du1[j] + dt3*du2[j] + dt3*du3[j] + dt6*du4[j];
	}

	// Advance time
	s->t += dt;
}





void Solver::init(Physical_parameters params, Control_parameters control, const std::vector<f_of_t_table_record>& l_table_functions_of_t,
	const std::vector<f_of_t_n_table_record>& l_table_functions_of_t_n)
{
	// Copy observables tables
	table_functions_of_t = l_table_functions_of_t;
	table_functions_of_t_n = l_table_functions_of_t_n;

	// Assign physical parameters
	N = params.N;
	F0 = params.F0;
	F1 = params.F1;
	nu = params.nu;

	// Assign control parameters
	dt = control.dt;
	ntimesteps = control.ntimesteps;
	calcsteps = control.calcsteps;
	dumpsteps = control.dumpsteps;

	// Reset observables and statistics
	current.init(N, table_functions_of_t, table_functions_of_t_n);
	current.to_zero();
	stat.init(N, table_functions_of_t, table_functions_of_t_n);
	stat.to_zero();
	stepcnt = 0;
	statcnt = 0;

	// Init solution
	s = new Solution(N);

    // Init arrays 
    k.resize    (N, 0.);
    k2.resize   (N, 0.);
	s->u.resize (N, C0);
    u2.resize   (N, C0);
	u3.resize   (N, C0);
	u4.resize   (N, C0);
    du1.resize  (N, C0);
    du2.resize  (N, C0);
	du3.resize  (N, C0);
	du4.resize  (N, C0);
    
    // Init wavenumbers
    for (int j=0; j<N; j++) {
		k[j] = K0 * pow(lam, j);
        k2[j] = k[j] * k[j]; 
    }
}

void Solver::assign_solution(Solution s0)
{
	check(s->u.size() == s0.u.size());
	for (unsigned int i = 0; i < s->u.size(); i++)
		s->u[i] = s0.u[i];
}

void Solver::assign_statistics(observables& stat0)
{
	check(stat0.get_N() == N);
	if (stat.get_N() != N)
	{
		stat.free();
		stat.init(N, table_functions_of_t, table_functions_of_t_n);
	}
	stat.to_zero();
	stat.sum(stat0);
}


void Solver::calculate_observables()
{// Observables calculation START
    
    if ((stepcnt % calcsteps)==0) {

        //double sum1=0, sum2=0, sum3=0, sum4=0, sum_k4u2=0;
        //double dd1=0, dd2=0, dd3=0;
        current.to_zero();
        
        
        for(int j=0; j<N; j++){
			double pm = (j % 2 == 0) ? (1.0) : (-1.0);
			double nuk2 = (nu)*k2[j];
			double pm_nuk3 = pm*nuk2*k[j];

			double u_sq = (double)real(s->u[j] * conj(s->u[j]));
			double u_abs = sqrt(u_sq);

            // Totals START
			current[Ot::Ot_Etot] += u_sq;
			current[Ot::Ot_Htot] += pm*k[j] * u_sq;
			current[Ot::Ot_Ediss] += -2.*nuk2*u_sq;
			current[Ot::Ot_Hdiss] += -2.*pm_nuk3*u_sq;
			current[Ot::Ot_Omega] += k2[j] * u_sq;
            if (j==0) {
				current[Ot::Ot_Ein] += 2.*real(s->u[0]*conj(F0));
				current[Ot::Ot_Hin] += 2.*pm*k[0]*real(s->u[0]*conj(F0));
            }
            if (j==1) {
				current[Ot::Ot_Ein] += 2.*real(s->u[1]*conj(F1));
				current[Ot::Ot_Hin] += 2.*pm*k[1]*real(s->u[1]*conj(F1));
            }

            // Totals END



			// Spectrals START
			current[Otn::Otn_REAL_Un][j] = real(s->u[j]);
			current[Otn::Otn_IMAG_Un][j] = imag(s->u[j]);
			current[Otn::Otn_En][j] = u_sq;

			// Fluxes START
			complex<double> C3n = C0;
			complex<double> C3np1 = C0;

			if (((j - 1) >= 0) && ((j + 1) < N))  C3n = s->u[j + 1] * conj(s->u[j] * s->u[j - 1]);
			if ((j + 2) < N)                      C3np1 = s->u[j + 2] * conj(s->u[j + 1] * s->u[j]);

			current[Otn::Otn_C3][j] = imag(C3n);
			current[Otn::Otn_Eflux][j] = (-2.*k[j]) * (A*lam*imag(C3np1) + (B + A)*imag(C3n));
			current[Otn::Otn_Hflux][j] = (-2.*pm*k[j] * k[j]) * (A*lam*imag(C3np1) + (B - (A / lam))*imag(C3n));
			// Fluxes END


			// Velocity structure functions START
			current[Otn::Otn_SF1][j] = u_abs;
			current[Otn::Otn_SF2][j] = u_sq;
			current[Otn::Otn_SF3][j] = current[Otn::Otn_SF2][j] * current[Otn::Otn_SF1][j];
			current[Otn::Otn_SF4][j] = current[Otn::Otn_SF3][j] * current[Otn::Otn_SF1][j];
			current[Otn::Otn_SF5][j] = current[Otn::Otn_SF4][j] * current[Otn::Otn_SF1][j];
			current[Otn::Otn_SF6][j] = current[Otn::Otn_SF5][j] * current[Otn::Otn_SF1][j];
			// Velocity structure functions END
			// Spectrals END
        }  // j
        

        // cumulate statistics
		stat.sum(current);
		++statcnt;
    }
    
}// Observables calculation END



Solution::Solution(unsigned int l_N)
{
	u.resize(l_N, complex<double>(0., 0.));
	t = 0.;
}

Solution::Solution(unsigned int l_N, unsigned long int random_seed)
{
	u.resize(l_N, complex<double>(0., 0.));
	
	srand(random_seed);

	unsigned int nmax = (default_initialization_n_max > l_N) ? (l_N) : (default_initialization_n_max);

	double k = initialization_K0;
	for (int j = 0; j < nmax; j++)
	{
		double dtmp = pow(k, -(1. / 3.));
		double ra = TWOPI*((double)rand() / (double)RAND_MAX);
		u[j] = dtmp*complex<double>(cos(ra), sin(ra));
		k *= initialization_lambda;
	}

	t = 0.;
}



Solution::Solution(unsigned int l_N, std::string filename)
{
	u.resize(l_N);

	// Check if input file (for velocity) exists
	File in(filename);
	in.open("r");

	// Load
	double v1r, v1i;
	cout << "# reading from input file: " << in.get_name() << "\n";
	for (int j = 0; j < l_N; j++) {
		if (feof(in.cfile())) {
			cout << "#WARNING! initial configuration has less shells than the system. Remaining shells will be padded to zero.\n";
			break;
		}
		int jj;

		int ee = fscanf(in.cfile(), "%d %le %le\n", &jj, &v1r, &v1i);
		check(ee == 3);
		u[j] = complex<double>(v1r, v1i);
		cout << j << " : " << u[j] << "\n";

	}
	if (!feof(in.cfile())) {
		cout << "#WARNING! initial configuration has more shells than the system. Shellnumbers >= "<<l_N<<" will be truncated.\n";
	}
	in.close();

	t = 0.;
}

Solution::~Solution()
{
	if (!u.empty()) u.clear();
}
