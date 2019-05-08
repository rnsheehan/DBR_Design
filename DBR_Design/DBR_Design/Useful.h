#ifndef USEFUL_H
#define USEFUL_H

// Library of functions that are very useful
// R. Sheehan 4 - 7 - 2011

namespace useful_funcs{
	
	std::string TheTime();

	void exit_failure_output(std::string reason);

	bool valid_filename_length(const std::string &name);
}

// This class will describe the bounds of a parameter sweep space
// e.g. sweep from start-value to stop-value in steps of delta; 
// R. Sheehan 1 - 3 - 2019

class sweep {
public:
	sweep(); 
	sweep(int &n_pts, double &start_val, double &stop_val);
	sweep(sweep &swp_obj); 
	~sweep(); 

	// setter
	void set_vals(int &n_pts, double &start_val, double &stop_val);
	void set_vals(sweep &swp_obj); 

	// getters
	inline int get_Nsteps() { return Nsteps;  }

	inline bool defined() { return params_defined;  }

	inline double get_start() { return start; }
	inline double get_stop() { return stop; }
	inline double get_delta() { return delta; }

	double get_val(int i); // access the ith computed sweep parameter

	std::vector<double> get_vals(); // get all the computed sweep parameters
private:
	int Nsteps; 

	bool params_defined; 

	double stop; 
	double start; 
	double delta; 

	std::vector<double> values; 
};

#endif