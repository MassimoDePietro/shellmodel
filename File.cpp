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


/// \file File.cpp
/// \brief This file contains the implementation of the File class.



#include "File.hpp"
using namespace std;


void File::open(const char* mode)
{
	if (m_is_open) return;
	if (m_fname == "") exit_with_message("ERROR! You tried to open a file without specifying the filename.");
	
	m_fp = fopen(m_fname.c_str(), mode); 
	if (m_fp == nullptr) exit_with_message(string("ERROR! Could not open file \"") + m_fname + string("\"\n"));

	m_is_open = true;
}


void File::close()
{
	if (!m_is_open) return;
	fclose(m_fp);
	m_fp = nullptr;
	m_is_open = false;
}

FILE* File::cfile() { return m_fp; }

string& File::get_name() { return m_fname; }

void File::set_name(const string& fname) { m_fname = fname; }

bool File::is_open()  {return m_is_open;}

