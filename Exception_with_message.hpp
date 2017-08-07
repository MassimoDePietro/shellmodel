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


/// \file Exception_with_message.hpp
/// \brief This file contains the exit_with_message() function, which prints an
/// error message and abort the program.


#pragma once

#include <cstdlib>
#include <string>

void exit_with_message(std::string message);
