#ifndef ATTACH_H
#include "Attach.h"
#endif

dispersion::dispersion()
{
	// Default Constructor
	params_defined = false;

	neff_calc = nullptr; core = nullptr; substrate = nullptr; cladding = nullptr; 
}

dispersion::dispersion(sweep &swp_obj, material *Ncore, material *Nsub, material *Nclad)
{
	set_params(swp_obj, Ncore, Nsub, Nclad);
}

dispersion::~dispersion()
{
	// Deconstructor

	clear(); 
}

void dispersion::set_params(sweep &swp_obj, material *Ncore, material *Nsub, material *Nclad)
{
	// Assign the values to the object parameters
	try {
		bool c1 = swp_obj.defined();
		bool c4 = Ncore != nullptr ? true : false; 
		bool c5 = Nsub != nullptr ? true : false; 
		bool c6 = Nclad != nullptr ? true : false; 
		bool c10 = c1 && c4 && c5 && c6; 

		if (c10) {

			wavelength.set_vals(swp_obj); // assign the values for the wavelength parameter space

			// assign the material objects
			core = Ncore; 
			substrate = Nsub; 
			cladding = Nclad; 

			neff_vals.clear(); 
			ng_vals.clear(); 

			params_defined = true; 
		}
		else {
			std::string reason; 
			reason = "Error: void dispersion::set_params(double wl_start, double wl_finish, double wl_step, material *Ncore, material *Nsub, material *Nclad)\n"; 
			if (!c1) reason += "swp_obj is not correct\n"; 
			if (!c4) reason += "Ncore has not been correctly assigned"; 
			if (!c5) reason += "Nsub has not been correctly assigned"; 
			if (!c6) reason += "Nclad has not been correctly assigned"; 
			throw std::invalid_argument(reason); 
		}

	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void dispersion::compute_dispersion(bool polarisation, wg_dims &dim_obj, std::string &filename, bool loud)
{
	// Compute the dispersion curve data based on the input parameters
	// R. Sheehan 7 - 3 - 2019

	try {
		bool c1 = params_defined ? true : false; 
		bool c2 = wavelength.defined() ? true : false; 
		bool c3 = dim_obj.defined() ? true : false; 
		bool c4 = neff_calc != nullptr ? true : false; 
		bool c10 = c1 && c2 && c3 && c4; 

		if (c10) {
			// Compute dispersion curves here

			neff_vals.clear();	ng_vals.clear();

			double lambda, ncore, nclad, nsub; // variable for storing wavelength value

			ri_vals ri_obj; // Declarate Refractive Index object

			for (int i = 0; i < wavelength.get_Nsteps(); i++) {
				lambda = wavelength.get_val(i); // Assign the current value of the wavelength

				core->set_wavelength(lambda); // tell the refractive index objects what the wavelength is

				substrate->set_wavelength(lambda);

				cladding->set_wavelength(lambda);

				// Assign the RI values to the RI object
				ncore = core->refractive_index(); 
				nsub = substrate->refractive_index(); 
				nclad = cladding->refractive_index(); 
				ri_obj.set_rib_wire(ncore, nsub, nclad, lambda);

				// Assign the calculation parameters to the EIM object
				neff_calc->set_params(polarisation, dim_obj, ri_obj);

				// perform the EIM calculation
				neff_calc->reduce_wg();

				neff_calc->get_index(false);

				neff_vals.push_back( neff_calc->neff_value() );

				if(loud) std::cout << lambda << " , " <<	ncore <<" , " << nsub << " , " << nclad << " , " <<	neff_calc->neff_value() << "\n"; 
			}

			compute_group_index(loud); // compute the group index data from the effective index data

			save_data_to_file(filename); // save neff and ngroup data to file
		}
		else {
			std::string reason;
			reason = "Error: void dispersion::compute_dispersion()\n";
			if(!c1) reason += "Parameters not defined\n"; 
			if(!c2) reason += "wavelength sweep params not defined\n"; 
			if(!c3) reason += "waveguide dimensions not defined\n"; 
			if(!c4) reason += "EIM calculation object not defined\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void dispersion::compute_group_index(bool loud)
{
	// Compute the group index data from a set of effective index data
	// R. Sheehan 15 - 4 - 2019

	try {
		bool c1 = wavelength.defined() ? true : false;
		bool c2 = static_cast<int>(neff_vals.size()) == wavelength.get_Nsteps() ? true : false; 
		bool c10 = c1 && c2; 

		if (c10) {
			ng_vals.clear(); 

			int N = static_cast<int>(neff_vals.size()); 
			double df, dx, ng_val; 
			for (int i = 0; i < N; i++) {
				if (i == 0) {
					// compute n'(lambda) using forward finite difference
					df = neff_vals[i + 1] - neff_vals[i]; 
					dx = wavelength.get_val(i + 1) - wavelength.get_val(i); 
				}
				else if (i == N - 1) {
					// compute n'(lambda) using backward finite difference
					df = neff_vals[i] - neff_vals[i-1];
					dx = wavelength.get_val(i) - wavelength.get_val(i-1);
				}
				else {
					// compute n'(lambda) using central finite difference
					df = neff_vals[i + 1] - neff_vals[i - 1];
					dx = wavelength.get_val(i + 1) - wavelength.get_val(i-1);
				}
				
				// compute group index from n_{g} = n_{eff} - lambda d n_{eff} / d lambda
				if (dx > 0) {
					ng_val = neff_vals[i] - wavelength.get_val(i)*(df / dx); 
				}
				else {
					ng_val = 0.0; 
				}

				ng_vals.push_back(ng_val); 

				if (loud) std::cout << wavelength.get_val(i)<<" , "<< neff_vals[i] << " , " << ng_val << "\n"; 
			}
		}
		else {
			std::string reason;
			reason = "Error: void dispersion::compute_group_index()\n";
			if (!c1) reason += "wavelength sweep params not defined\n";
			if (!c2) reason += "neff data not defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void dispersion::save_data_to_file(std::string &filename)
{
	// Write the computed dispersion data to a file
	// R. Sheehan 15 - 4 - 2019

	try {
		bool c1 = wavelength.defined() ? true : false;
		bool c2 = static_cast<int>(neff_vals.size()) == wavelength.get_Nsteps() ? true : false;
		bool c3 = static_cast<int>(ng_vals.size()) == wavelength.get_Nsteps() ? true : false;
		bool c4 = filename != empty_str ? true : false; 
		bool c5 = useful_funcs::valid_filename_length(filename); 
		bool c10 = c1 && c2 && c3 && c4 && c5;

		if (c10) {
			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

			if (write.is_open()) {
				for (int i = 0; i < wavelength.get_Nsteps(); i++) {
					write << std::setprecision(10) << wavelength.get_val(i) << " , " << neff_vals[i] << " , " << ng_vals[i] << "\n"; 
				}

				write.close(); 
			}
			else {
				std::string reason;
				reason = "Error: void dispersion::save_data_to_file(std::string &filename)\n";
				reason += "Could not open file: " + filename + "\n";
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void dispersion::save_data_to_file()\n";
			if (!c1) reason += "wavelength sweep params not defined\n";
			if (!c2) reason += "neff data not defined\n";
			if (!c3) reason += "ng data not defined\n";
			if (!c4) reason += "filename not defined\n";
			if (!c5) reason += "filename: " + filename + " too long\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what();
	}
}

std::vector<std::vector<double>> dispersion::get_computed_data(bool cnvrt_wl_nm)
{
	// return computed dispersion data as the rows of a new array
	// R. Sheehan 8 - 5 - 2019

	try {
		bool c1 = wavelength.defined() ? true : false;
		bool c2 = static_cast<int>(neff_vals.size()) == wavelength.get_Nsteps() ? true : false;
		bool c3 = static_cast<int>(ng_vals.size()) == wavelength.get_Nsteps() ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::vector<std::vector<double>> data; 

			data.push_back(wavelength.get_vals()); 
			data.push_back(neff_vals); 
			data.push_back(ng_vals); 

			if (cnvrt_wl_nm) {
				// convert wavelength data to nm
				for (size_t i = 0; i < data[0].size(); i++) {
					data[0][i] *= 1000.0; 
				}
			}

			return data; 
		}
		else {
			return std::vector<std::vector<double>>(); 
			std::string reason;
			reason = "Error: std::vector<std::vector<double>> dispersion::get_computed_data()\n";
			if (!c1) reason += "wavelength sweep params not defined\n";
			if (!c2) reason += "neff data not defined\n";
			if (!c3) reason += "ng data not defined\n";
			throw std::invalid_argument(reason);
		}	
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void dispersion::clear()
{
	params_defined = false;

	neff_vals.clear();	ng_vals.clear();
}

// Defintion of the wire eaveguide dispersion calculation
wire_dispersion::wire_dispersion()
{
	// Default constructor
}

void wire_dispersion::compute_dispersion_data(bool polarisation, sweep &swp_obj, wg_dims &dim_obj, material *Ncore, material *Nsub, material *Nclad, std::string &filename, bool loud)
{
	// Compute the dispersion data based on the defined inputs
	// test to ensure that input dim_obj is associated with a wire waveguide

	try {
		if (dim_obj.get_wg_code() == WIRE_WG) {	 
			
			Wire w_obj; // Declarate Wire waveguide object

			neff_calc = &w_obj; // Point the EIM object to the Wire waveguide object

			set_params(swp_obj, Ncore, Nsub, Nclad);				

			compute_dispersion(polarisation, dim_obj, filename, loud); 
		}
		else {
			std::string reason;
			reason = "Error: void wire_dispersion::compute_dispersion_data()\n";
			reason += "Input dimension values are not correct\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Definition of the rib waveguide dispersion calculation

rib_dispersion::rib_dispersion()
{
	// Default Constructor
}

void rib_dispersion::compute_dispersion_data(bool polarisation, sweep &swp_obj, wg_dims &dim_obj, material *Ncore, material *Nsub, material *Nclad, std::string &filename, bool loud)
{
	// Compute the dispersion data based on the defined inputs
	// test to ensure that input dim_obj is associated with a wire waveguide

	try {
		if (dim_obj.get_wg_code() == RIB_WG) {
			// Compute dispersion curves here
			Rib w_obj; // Declarate Rib waveguide object

			neff_calc = &w_obj; // Point the EIM object to the Wire waveguide object

			set_params(swp_obj, Ncore, Nsub, Nclad);

			compute_dispersion(polarisation, dim_obj, filename, loud);
		}
		else {
			std::string reason;
			reason = "Error: void rib_dispersion::compute_dispersion_data()\n";
			reason += "Input dimension values are not correct\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}