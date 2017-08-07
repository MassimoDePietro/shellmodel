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


/// \file main.cpp
/// \brief All the stuff pertaining the interaction with the system is in this file.
/// Function main() is here. 



#include <iostream>
#include <map>

//#include "Exception_with_message.hpp"
#include "check.hpp"
//#include "Command_line_arguments.hpp"
#include "File.hpp"
#include "Observables.hpp"
#include "Solver.hpp"



using namespace std;



/// A simple class for storing command line arguments.
struct Command_line_argument {
	Command_line_argument(string n, string d, string v) : 
		name{ n }, 
		description{ d }, 
		value{ v } 
	{};

	string name;			///< Name of the argument.
	string description;		///< Description of the argument.
	string value;			///< Value of the argument (kept as string).
};

/// A simple class for storing options regarding the execution of the program
/// and its interaction with the system.
struct Program_options {
	bool flag_load_u0 = false;			///< Load initial velocity from file
										///< (default: false)
	bool flag_load_statistics = false;	///< Load initial statistics from file
										///< (default: false)

	std::string input_prefix{ "start" };	///< Name of the input file
	std::string output_prefix{ "sim" };		///< Name of the output file

	long long random_seed{ 123456 };	///< Seed for the Random Number Generator
};

/// A simple class collecting the files used for input.
struct Input_files {
	File F_in_u   = File("start");		///< Name of the input file for velocity
	File F_in_obs = File("start-obs");	///< Name of the input file for statistics
};

/// A simple class collecting the files used for output.
struct Output_files {
	File F_out_Ot  = File("sim-totals");	///< Name of the output file for obsrvables of type f_of_t
	File F_out_Otn = File("sim-spectrals"); ///< Name of the output file for obsrvables of type f_of_t_n
	File F_out_u   = File("sim-last");		///< Name of the output file for the last velocity configuration
	File F_out_obs = File("sim-last-obs");	///< Name of the output file for the statistics
};


//----------------------------------------//
// ---------- Global Functions ---------- //
//----------------------------------------//

/// This function is called at the end of the simulation to write the file used 
/// for restart.
void dump_last(Solver& s, Output_files& output_files)
{
    // file last
    File& fLast = output_files.F_out_u;
    fLast.open("w");
    for(int j=0; j<s.N; j++){
        fprintf(fLast.cfile(), "%d %1.16e %1.16e\n",j, 
			real(s.s->u[j]), imag(s.s->u[j]));
    }
    fLast.close();
    
    
    // file last-observables
    File& fLastObs = output_files.F_out_obs;
    fLastObs.open("wb");
    fwrite(&(s.s->t), sizeof(s.s->t), 1, fLastObs.cfile());
    fwrite(&(s.statcnt), sizeof(s.statcnt), 1, fLastObs.cfile());
	s.current.save_to_stream_binary(fLastObs.cfile());
	s.stat.save_to_stream_binary(fLastObs.cfile());
    fLastObs.close();
}



/// This function is called during the execution to dump on files the values of 
/// the observables (both instantaneous and running averages).
void dump_during_run(Solver& s, Output_files& output_files)
{
    FILE* ftotals = output_files.F_out_Ot.cfile();
    FILE* fspectrals = output_files.F_out_Otn.cfile();

    if ((s.stepcnt % s.dumpsteps) == 0) {
        double sden = 1./(double)s.statcnt;

        // Dump total observables (summed on all shells)
        fprintf(ftotals, "%e ", s.s->t);							   // 1 time
        for(int i=0; i < s.table_functions_of_t.size(); i++) {
			auto obs_id = s.table_functions_of_t[i].type;
            fprintf(ftotals, "%e ", s.current[obs_id]);
            fprintf(ftotals, "%e ", s.stat[obs_id]*sden);
        }
        fprintf(ftotals, "\n");
        fflush(ftotals);
        
		
        // Dump spectral observables
        for(int j=0; j<s.N; j++) {
            fprintf(fspectrals, "%e ", s.s->t);                        // 1 time
            fprintf(fspectrals, "%d ", j);                             // 2 n
			for (int i = 0; i < s.table_functions_of_t_n.size(); i++) {
				auto obs_id = s.table_functions_of_t_n[i].type;
				fprintf(fspectrals, "%e ", (s.current[obs_id])[j]);
				fprintf(fspectrals, "%e ", (s.stat[obs_id])[j] * sden);
			}
            fprintf(fspectrals, "\n");
        }
		fprintf(fspectrals, "\n\n");
        fflush(fspectrals);
    }
}





/// This is one of the initialization functions. It prepares the files and makes
/// sure there are no problems opening/closing them.
void prepare_files(Program_options& options, Input_files& input_files, 
	Output_files& output_files, 
	const vector<f_of_t_table_record>& table_functions_of_t, 
	const vector<f_of_t_n_table_record>& table_functions_of_t_n)
{
	// Fill in filenames
	 input_files.F_in_u.set_name(options.input_prefix);
	 input_files.F_in_obs.set_name(options.input_prefix + "-obs");
	output_files.F_out_Ot.set_name(options.output_prefix + "-totals");
	output_files.F_out_Otn.set_name(options.output_prefix + "-spectrals");
	output_files.F_out_u.set_name(options.output_prefix + "-last");
	output_files.F_out_obs.set_name(options.output_prefix + "-last-obs");

	// Try opening all output files, check for errors
	vector<File> flist({ output_files.F_out_Ot, output_files.F_out_Otn, 
		                 output_files.F_out_u, output_files.F_out_obs });
	for (auto& f : flist) {
		f.open("w");
		f.close();
	}

	// Keep open the totals and spectrals files
	output_files.F_out_Ot.open("w");
	output_files.F_out_Otn.open("w");

	// Write columns descriptor on the totals and spectrals files
	FILE* ff;
	int col;

	// file -totals
	ff = output_files.F_out_Ot.cfile();
	fprintf(ff, "# %d : %s\n", 1, "time");
	col = 2;
	//for(int i=0; i < Ot::MAX; i++) {
	for (auto& obs : table_functions_of_t) {
		const auto& name = obs.name.c_str();
		fprintf(ff, "# %d : %s\n", col, name);
		fprintf(ff, "# %d : %s running avg\n", col + 1, name);
		col += 2;
	}

	// file -spectrals
	ff = output_files.F_out_Otn.cfile();
	fprintf(ff, "# %d : %s\n", 1, "time");
	fprintf(ff, "# %d : %s\n", 2, "n");
	col = 3;
	for (auto& obs : table_functions_of_t_n) {
		const auto& name = obs.name.c_str();
		fprintf(ff, "# %d : %s\n", col, name);
		fprintf(ff, "# %d : %s running avg\n", col + 1, name);
		col += 2;
	}
}









/// This is one of the initialization functions. It reads the command line 
/// arguments into a vector of Command_line_argument.
void parse_args(int argc, char *argv[], vector<Command_line_argument>& cla)
{
	// Register the command line arguments needed
	cla.emplace_back("<load>", "0 for starting new simulation; 1 for loading only the initial velocity configuration; 2 for loading both an initial velocity configuration and statistics.", string());
	cla.emplace_back("<input_file_prefix>", "Name of the input file for the initial configuration. The file with the initial statistics, if used, must be named <input_file_prefix>-obs.", string());
	cla.emplace_back("<output_file_prefix>", "The simulation will produce the output files: <output_file_prefix>-totals, <output_file_prefix>-spectrals, <output_file_prefix>-last, <output_file_prefix>-last-obs.", string());
	cla.emplace_back("<N>", "Number of shells.", string());
	cla.emplace_back("<absF0>", "Magnitude of forcing (on shell 0 only).", string());
	cla.emplace_back("<argF0>", "Phase of forcing (on shell 0 only).", string());
	cla.emplace_back("<nu>", "Viscosity.", string());
	cla.emplace_back("<dt>", "Integration timestep.", string());
	cla.emplace_back("<ntimesteps>", "Number of timesteps to integrate", string());
	cla.emplace_back("<calcsteps>", "Observables will be calculated every <calcsteps> timesteps. Statistics will be updated with the same frequency.", string());
	cla.emplace_back("<dumpsteps>", "Observables will be dumped every <dumpsteps> timesteps.", string());
	cla.emplace_back("<random_seed>", "Seed for the random number generator. Currently used only when generating an initial condition.", string());
	

	// Check number of arguments, eventually print usage
	if ((argc - 1) != cla.size()) {
		cerr << "Usage:\n";
		for (auto& arg : cla) cout << arg.name << " ";
		cout << "\n\n";
		cout << "Detailed explaination of all the arguments:\n";
		for (auto& arg : cla) {
			cout << "\t" << arg.name << ": " << arg.description << "\n";
		}
	}

	// Read the command line arguments
	for (unsigned int i = 0; i < cla.size(); i++) cla[i].value = argv[i+1];
}



/// This is one of the initialization functions. It parses the read 
/// Command_line_argument's and put the values of the parameters in the appropriate containers.
void fill_in_parameters(vector<Command_line_argument>& cla, Program_options& options, Physical_parameters& params, Control_parameters& control)
{
	// Create a map for convenience
	map<string, string> args;
	for (auto& arg : cla) args.emplace(arg.name, arg.value);

	// Fill in options
	int pload = stoi(args["<load>"]);
	if ((pload < 0) || (pload > 2)) exit_with_message("<load> must be 0,1 or 2");
	if (pload > 0)  options.flag_load_u0 = true;
	if (pload == 2) options.flag_load_statistics = true;
	options.random_seed = stoll(args["<random_seed>"]);
	options.input_prefix = args["<input_file_prefix>"];
	options.output_prefix = args["<output_file_prefix>"];

	// Fill in physical parameters
	params.N = static_cast<unsigned int>(stoi(args["<N>"]));
	if (params.N < 3) exit_with_message("It must be N >= 3");

	double rho = stod(args["<absF0>"]);
	double alpha = stod(args["<argF0>"]);
	if (rho < 0.0) exit_with_message("It must be <absF0> >= 0");
	params.F0 = complex<double>(rho*cos(alpha), rho*sin(alpha));
	params.F1 = complex<double>(0., 0.);
	params.nu = stod(args["<nu>"]);
	if (params.nu < 0.0) exit_with_message("It must be <nu> >= 0");

	// Fill in integration and control parameters
	control.dt = stod(args["<dt>"]);
	control.ntimesteps = stoull(args["<ntimesteps>"]);
	control.calcsteps = stoull(args["<calcsteps>"]);
	control.dumpsteps = stoull(args["<dumpsteps>"]);
	if (control.ntimesteps < 0) exit_with_message("It must be <ntimesteps> >= 0");
	if (control.calcsteps < 0) exit_with_message("It must be <calcsteps> >= 0");
	if (control.dumpsteps < 0) exit_with_message("It must be <dumpsteps> >= 0");
	if (control.ntimesteps < control.calcsteps) exit_with_message("It must be <ntimesteps> >= <calcsteps>");
	if (control.ntimesteps < control.dumpsteps) exit_with_message("It must be <ntimesteps> >= <dumpsteps>");
	if (control.dumpsteps < control.calcsteps) exit_with_message("It must be <dumpsteps> >= <calcsteps>");
}



/// This function prints on screen all the values of the parameters of the system,
/// as a check.
void print_options_and_parameters(Program_options& options, Physical_parameters& params, Control_parameters& control, Input_files& input_files, Output_files& output_files)
{
	// Options
	cout << "# Options:\n";
	cout << "# -- flag_load_u0 = " << options.flag_load_u0 << "\n";
	cout << "# -- flag_load_statistics = " << options.flag_load_statistics << "\n";
	cout << "# -- input_prefix = " << options.input_prefix << "\n";
	cout << "# -- output_prefix = " << options.output_prefix << "\n";
	cout << "# -- random_seed = " << options.random_seed << "\n";

	if (options.flag_load_u0) cout << "# Input (velocity) will be read from " << input_files.F_in_u.get_name() << "\n";
	if (options.flag_load_statistics) cout << "# Input (statistics) will be read from " << input_files.F_in_obs.get_name() << "\n";

	cout << "# Output will be written on files:\n";
	vector<File> flist({ output_files.F_out_Ot, output_files.F_out_Otn, output_files.F_out_u, output_files.F_out_obs });
	for (auto& of : flist) {
		cout << "# " << of.get_name() << "\n";
	}

	// Physical parameters
	cout << "# Physical parameters:\n";
	cout << "# -- N = " << params.N << "\n";
	cout << "# -- F0 = " << params.F0 << "\n";
	cout << "# -- F1 = " << params.F1 << "\n";
	cout << "# -- nu = " << params.nu << "\n";
	cout << "# -- (hardcoded) Lambda = " << params.lam << "\n";
	cout << "# -- (hardcoded) K0 = " << params.K0 << "\n";
	cout << "# -- (hardcoded) A = " << params.A << "\n";
	cout << "# -- (hardcoded) B = " << params.B << "\n";
	cout << "# -- (hardcoded) C = " << params.C << "\n";

	// Control
	cout << "# Control parameters:\n";
	cout << "# -- dt = " << control.dt << "\n";
	cout << "# -- ntimesteps = " << control.ntimesteps << "\n";
	cout << "# -- calcsteps = " << control.calcsteps << "\n";
	cout << "# -- dumpsteps = " << control.dumpsteps << "\n";
	cout << "# Total integration time = " << control.dt*(double)control.ntimesteps << "\n";
	cout << "# Calculations every = " << control.dt*(double)control.calcsteps << "\n";
	cout << "# Statistics will contain # elements = " << control.ntimesteps / control.calcsteps << "\n";
	cout << "# Dump every = " << control.dt*(double)control.dumpsteps << "\n";
	cout << "# There will be # dumps = " << control.ntimesteps / control.dumpsteps << "\n";
}








/// \brief This in the main() function. Prepares and directs the shell model 
///  solver and its interaction with the system.
///
/// The structure is the following: 
/// - definition of the tables of the observables used;
/// - allocation of all the containers of the parameters etc.;
/// - parse command line arguments and store them in the containers of the parameters;
/// - alloc and initialize the shell model solver with the parameters;
/// - if necessary, read initial velocity and statistics and pass them to the shell model solver;
/// - initialization finished: dump initial values of the observables;
/// - MAIN LOOP:
///   + advance one timestep;
///   + calculate observables, update statistics, dump output;
/// - simulation finished: write files for restart;
/// - clean and exit.
///
int main(int argc, char *argv[])
{
	// Define observables and files to be used
	const vector<f_of_t_table_record> table_functions_of_t
	{
		// Table with columns:
		// Ot, name (as it will appear in the output files)
		{ Ot::Ot_Etot,  "(Total energy)" },
		{ Ot::Ot_Htot,  "(Total helicity)" },
		{ Ot::Ot_Ediss, "(Total energy dissipation)" },
		{ Ot::Ot_Hdiss, "(Total helicity dissipation)" },
		{ Ot::Ot_Ein,   "(Total energy injection)" },
		{ Ot::Ot_Hin,   "(Total helicity injection)" },
		{ Ot::Ot_Omega, "(Total enstrophy)" }
	};

	const vector<f_of_t_n_table_record> table_functions_of_t_n
	{
		// Table with columns:
		// Ot, name (as it will appear in the output files)
		{ Otn::Otn_REAL_Un, "(Re[u_n])" },
		{ Otn::Otn_IMAG_Un, "(Im[u_n])" },
		{ Otn::Otn_En,      "(Energy spectrum)" },
		{ Otn::Otn_Eflux,   "(Energy flux)" },
		{ Otn::Otn_Hflux,   "(Helicity flux)" },
		{ Otn::Otn_C3,      "(Third order correlation function)" },
		{ Otn::Otn_SF1,     "(Structure function, order 1)" },
		{ Otn::Otn_SF2,     "(Structure function, order 2)" },
		{ Otn::Otn_SF3,     "(Structure function, order 3)" },
		{ Otn::Otn_SF4,     "(Structure function, order 4)" },
		{ Otn::Otn_SF5,     "(Structure function, order 5)" },
		{ Otn::Otn_SF6,     "(Structure function, order 6)" }
	};


	vector<Command_line_argument> cla;
	Program_options options;
	Physical_parameters params;
	Control_parameters control;
	Input_files input_files;
	Output_files output_files;

	// Read command line arguments
	parse_args(argc, argv, cla);
	
	// Check and fill in parameters and options
	fill_in_parameters(cla, options, params, control);

	// Prepare files
	prepare_files(options, input_files, output_files, table_functions_of_t, table_functions_of_t_n);

    // Init solver
	Solver solver;
	solver.init(params, control, table_functions_of_t, table_functions_of_t_n);

	// Initial condition for the solution
	if (options.flag_load_u0) 
		solver.assign_solution(Solution(params.N, input_files.F_in_u.get_name()));
	else
		solver.assign_solution(Solution(params.N, options.random_seed));

	// Initial statistics
	if (options.flag_load_statistics) {
		File& f = input_files.F_in_obs;
		f.open("rb");
		cout << "# reading observables from input file: " << f.get_name() << "\n";
		observables o0;
		o0.load_from_stream_binary(f.cfile());
		f.close();
		solver.assign_statistics(o0);
	}

	// Print options and parameters
	print_options_and_parameters(options, params, control, input_files, output_files);

	// dump initial configuration on screen
	cout << "Initial config:\n";
	for (int j = 0; j < solver.N; j++) cout << j << " :" << solver.s->u[j] << "\n";

	// Write observables at t = t_initial
	solver.nlt();						// Calculate rhs for the observables calculations on the first timestep
	solver.calculate_observables();
	dump_during_run(solver, output_files);

    // MAIN LOOP
    for (long long stat1=1; stat1<= solver.ntimesteps; stat1++)
    {   
        // time step
		solver.advance_one_step_rk4();
    
		solver.calculate_observables();
        dump_during_run(solver, output_files);
    }
    
    // dump last configuration for restart
    dump_last(solver, output_files);

    // dump last configuration on screen
    cout << "Last config:\n";
	for (int j = 0; j < solver.N; j++) cout << j << " :" << solver.s->u[j] << "\n";
    
    // Clean and exit
	output_files.F_out_Ot.close();
	output_files.F_out_Otn.close();
    exit(0);
}



