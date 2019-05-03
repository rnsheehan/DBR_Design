#ifndef ATTACH_H
#include "Attach.h"
#endif

std::string useful_funcs::TheTime()
{
	// Implementation of a function returning the current time as a string
	// This is just a neater way of ensuring that the time can be correctly and easily accessed
	// without being hassled about whether or not you've remembered to use non-deprecated versions 
	// of certain functions
	// R. Sheehan 4 - 7 - 2011
	
	const int N=30;	
	char time_str[N];	
	size_t bytes=( N*sizeof(char) );
	
	time_t rawtime;
	
	struct tm timeinfo;
	struct tm *timeinfo_ptr;
	
	timeinfo_ptr=&timeinfo;
	
	// Get current time information
	time(&rawtime);
	
	localtime_s(timeinfo_ptr,&rawtime);
	
	asctime_s(time_str,bytes,timeinfo_ptr);
	
	// Deprecated calls
	//timeinfo=localtime(&rawtime);
	//asctime(timeinfo);
	
	std::string the_time;
	the_time.append(time_str);
	
	return the_time;
}

void useful_funcs::exit_failure_output(std::string reason)
{
	// Code that creates a file and writes a reason in it why the program crashed
	// If it is called of course
	// Call before using the exit(EXIT_FAILURE) command

	// This function outputs to a file an explanation of why the program exited with an EXIT_FAILURE
	// R. Sheehan 17 - 5 - 2011
	
	// Get current time information
	std::string time = TheTime();

	std::ofstream write; // open file for writing
	
	write.open("Exit_Failure_Explanation.txt",std::ios_base::out|std::ios_base::trunc);
	
	//if(!write){
	//	std::cout<<"You're not going to see this statement\n";
	//	std::cout<<"\n";
	//}
	//else{
	//	//printf ( "Current local time and date: %s", asctime (timeinfo) );
	//	write<<"Program Exit Explanation\n\n";
	//	write<<"Error occurred "<<time<<"\n";
	//	write<<reason<<"\n";
	//	write.close();
	//}

	if( write.is_open() ){
		
		write<<"Program Exit Explanation\n\n";
		write<<"Error occurred: "<<time<<"\n";
		write<<reason<<"\n";

		write.close();
	}
}

bool useful_funcs::valid_filename_length(const std::string &name)
{
	// Check that a string length is less than the MAX_PATH_LENGTH
	// This only really applies to windows

	return static_cast<int>(name.length()) < MAX_PATH_LENGTH ? true : false;
}

// Definition of the parameter sweep space class
// R. Sheehan 1 - 3 - 2019

sweep::sweep()
{
	// Default Constructor
	params_defined = false; 

	Nsteps = 0; 

	start = stop = delta = 0.0; 
}

sweep::sweep(int &n_pts, double &start_val, double &stop_val)
{
	// Primary Constructor

	set_vals(n_pts, start_val, stop_val); 
}

sweep::sweep(sweep &swp_obj)
{
	// Copy constructor
	set_vals(swp_obj); 
}

sweep::~sweep()
{
	// Deconstructor

	values.clear(); 
}

// setter
void sweep::set_vals(int &n_pts, double &start_val, double &stop_val)
{
	// Define parameters for the sweep object
	// R. Sheehan 1 - 3 - 2019

	// On the differences between the different cast methods in C++
	// https://stackoverflow.com/questions/332030/when-should-static-cast-dynamic-cast-const-cast-and-reinterpret-cast-be-used
	// Main message: Use static_cast for ordinary type conversions

	try {
		bool c1 = n_pts > 3 ? true : false; 
		bool c2 = fabs(stop_val - start_val) > 0 ? true : false;
		bool c10 = c1 && c2; 

		if (c10) {
			// Assign values to the parameters
			Nsteps = n_pts; 
			start = std::min(stop_val, start_val); 
			stop = std::max(stop_val, start_val); 
			delta = (stop - start) / (static_cast<double>(Nsteps - 1)); 

			// Fill the vector with the desired values
			values.clear(); 

			double pos = start; 
			for (int i = 0; i < Nsteps; i++) {
				values.push_back(pos);
				pos += delta; 
			}

			params_defined  = true; 
		}
		else {
			std::string reason;
			reason = "Error: void sweep::set_vals(int &n_pts, double &start_val, double &stop_val)\n";
			if (!c1) reason += "n_pts: " + template_funcs::toString(n_pts) + " is not correct\n";
			if (!c2) reason += "fabs(stop_val - start_val): " + template_funcs::toString(fabs(stop_val - start_val), 2) + " is not correct\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void sweep::set_vals(sweep &swp_obj)
{
	// Assign the values from one swp_obj to another

	try {
		if (swp_obj.defined()) {
			*this = swp_obj; 
		}
		else {
			std::string reason;
			reason = "Error: void sweep::void set_vals(sweep &swp_obj)\n";
			reason += "Parameters not defined for input object\n"; 

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double sweep::get_val(int i)
{
	// Access the computed sweep parameters stored in the array

	try {
		bool c1 = i > -1 ? true : false; 
		bool c2 = i < Nsteps ? true : false; 
		bool c3 = static_cast<int>(values.size()) == Nsteps ? true : false; 
		bool c4 = static_cast<int>(values.size()) > 3 ? true : false; 
		bool c10 = c1 && c2 && c3 && c4;

		if (c10) {
			return values[i]; 
		}
		else {
			return 0.0; 
			std::string reason; 
			reason = "Error: double sweep::get_val(int i)\n"; 
			if (!c1 || !c2) reason += "Out of bounds array access attempt\n"; 
			if (!c3 || !c4) reason += "values not correctly defined\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}
