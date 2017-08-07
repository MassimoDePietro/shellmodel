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


/// \file Observables.cpp
/// \brief This file contains the definition of the Observables class.


#include "observables.hpp"
#include "check.hpp"

#include <cassert>



using namespace std;


void observables::init(int N, const std::vector<f_of_t_table_record>& l_table_functions_of_t,
	const std::vector<f_of_t_n_table_record>& l_table_functions_of_t_n)
{
	this->N = N;

	// Copy tables
	table_functions_of_t = l_table_functions_of_t;
	table_functions_of_t_n = l_table_functions_of_t_n;

	// Calculate total size to allocate
	for (size_t el = 0; el < table_functions_of_t.size(); el++) num_f_of_t += 1;
	for (size_t el = 0; el < table_functions_of_t_n.size(); el++) num_f_of_t_n += N;

	f_of_t.resize(num_f_of_t, 0);
	f_of_t_n.resize(num_f_of_t_n, 0);
}

void observables::free()
{
	this->N = 0;
	this->num_f_of_t = 0;
	this->num_f_of_t_n = 0;
	
	f_of_t.clear();
	f_of_t_n.clear();
}

void observables::to_zero()
{
	for (auto& i : f_of_t) i = 0.;
	for (auto& i : f_of_t_n) i = 0.;
}

void observables::sum(const observables & obs)
{
	for (auto & i : table_functions_of_t) {
		auto observable_id = i.type;
		(*this)[observable_id] += obs[observable_id];
	}
	for (auto & i : table_functions_of_t_n) {
		auto observable_id = i.type;
		for (size_t j = 0; j < N; j++) {
			(*this)[observable_id][j] += obs[observable_id][j];
		}
	}
}

void observables::save_to_stream_binary(FILE * stream) const
{
	check(stream != nullptr);
	size_t first = std::fwrite(&N, sizeof(N), 1, stream);
	check(first == 1);
	size_t second = 0;
	for (auto& i : f_of_t) second += std::fwrite(&i, sizeof(double), 1, stream);
	for (auto& i : f_of_t_n) second += std::fwrite(&i, sizeof(double), 1, stream);
	check(second == num_f_of_t + num_f_of_t_n);
}

void observables::load_from_stream_binary(FILE * stream)
{
	check(stream != nullptr);
	size_t first = std::fread(&N, sizeof(N), 1, stream);
	check(first == 1);
	size_t second = 0;
	for (auto& i : f_of_t) second += std::fread(&i, sizeof(double), 1, stream);
	for (auto& i : f_of_t_n) second += std::fread(&i, sizeof(double), 1, stream);
	check(second == num_f_of_t + num_f_of_t_n);
}

unsigned int observables::get_N()
{
	return N;
}


double & observables::operator[](Ot x)
{
	return (f_of_t[static_cast<size_t>(x)]);
}

double* observables::operator[](Otn x)
{
	size_t address = static_cast<size_t>(x) * N;
	assert(address < f_of_t_n.size());
	return (&(f_of_t_n.data()[address]));
}

double observables::operator[](Ot x) const
{
	return (f_of_t[static_cast<size_t>(x)]);
}

const double * observables::operator[](Otn x) const
{
	size_t address = static_cast<size_t>(x) * N;
	assert(address < f_of_t_n.size());
	return (&(f_of_t_n.data()[address]));
}


void observables::add_to_statistics(const observables& current, observables& cumulative)
{
	cumulative.sum(current);
}

void observables::save_binary(FILE * stream, const observables & obs)
{
	obs.save_to_stream_binary(stream);
}

void observables::load_binary(FILE * stream, observables & obs)
{
	obs.load_from_stream_binary(stream);
}
