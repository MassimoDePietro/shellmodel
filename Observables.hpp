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


/// \file Observables.hpp
/// \brief This file contains the declaration of the Observables class and the 
/// enumerations defining all the observables used in the simulation. The 
/// Observables class stores and allows easy acces to all the observables 
/// measured in the simulation.

#pragma once

#include <vector>
#include <string>
#include <cstdio>




// Types of observables f_of_t
enum Ot : unsigned int {
	Ot_Etot = 0,
	Ot_Htot,
	Ot_Ediss,
	Ot_Hdiss,
	Ot_Ein,
	Ot_Hin,
	Ot_Omega
};
struct f_of_t_table_record {
	Ot type;
	std::string name;
};

// Types of observables f_of_t_n
enum Otn : unsigned int {
	Otn_REAL_Un = 0,
	Otn_IMAG_Un,
	Otn_En,
	Otn_Eflux,
	Otn_Hflux,
	Otn_C3,
	Otn_SF1,
	Otn_SF2,
	Otn_SF3,
	Otn_SF4,
	Otn_SF5,
	Otn_SF6
};
struct f_of_t_n_table_record {
	Otn type;
	std::string name;
};







class observables {
private:
	size_t N;				///< Number of shells
	size_t num_f_of_t;		///< How many f_of_t are there
	size_t num_f_of_t_n;	///< How many f_of_t_n are there
	std::vector<double> f_of_t;
	std::vector<double> f_of_t_n;

	std::vector<f_of_t_table_record> table_functions_of_t;
	std::vector<f_of_t_n_table_record> table_functions_of_t_n;

public:
	// Cosntructor/destructor
	observables() :
		N{ 0 },
		num_f_of_t{ 0 },
		num_f_of_t_n{ 0 }
	{};
	~observables() { this->free(); };

	// Methods
	void init(int N, const std::vector<f_of_t_table_record>& l_table_functions_of_t,
		const std::vector<f_of_t_n_table_record>& l_table_functions_of_t_n);
	void free();

	void to_zero();
	void sum(const observables& obs);

	void save_to_stream_binary(FILE* stream) const;
	void load_from_stream_binary(FILE* stream);

	unsigned int get_N();

	// Overload [] operator for an easier access to members
	double operator[] (Ot x) const;
	const double* operator[] (Otn x) const;
	double& operator[] (Ot x);
	double* operator[] (Otn x);


	// Static methods
	static void add_to_statistics(const observables& current, observables& cumulative);

	static void save_binary(FILE* stream, const observables& obs);
	static void load_binary(FILE* stream, observables& obs);
};

