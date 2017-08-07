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


/// \file File.hpp
/// \brief This file contains the declaration of the File class. The File class 
/// is a wrapper around the c FILE type with the addition that it stores the 
/// filename and checks for errors when opening.


#pragma once

#include <cstdio>
#include <string>
#include <vector>
#include <tuple>
#include "Exception_with_message.hpp"



class File {
public:
	File() {};
	File(const std::string& fname) : m_fname(fname) {};
	~File() {};
	
	void open(const char* mode);
	void close();
	
	FILE* cfile();
	std::string& get_name();
	void set_name(const std::string& fname);
	bool is_open();
	
private:
	FILE* m_fp = nullptr;
	bool m_is_open = false;
	std::string m_fname = "";
};



