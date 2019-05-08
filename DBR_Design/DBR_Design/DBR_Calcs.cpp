#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the methods declared in the DBR class
// R. Sheehan 7 - 5 - 2019

DBR::DBR()
{
	// Default constructor
	inputs_defined = false; coefficients_defined = false;
	n_wl_vals = n_periods = 0;
	Bragg_WL = Bragg_neff = Bragg_beta = Bragg_Period = min_feat_size = 0.0;
}

DBR::DBR(double &min_feature, double &wl_Bragg, std::vector<double> &wl_vals, std::vector<double> &neff_w_low, std::vector<double> &neff_w_mid, std::vector<double> &neff_w_high, std::vector<double> &ngroup_w_mid)
{
	// Constructor
	set_params(min_feature, wl_Bragg, wl_vals, neff_w_low, neff_w_mid, neff_w_high, neff_w_mid); 
}

DBR::~DBR()
{
	// Deconstructor
	wavelengths.clear(); 
	neff_m.clear(); neff_l.clear(); neff_h.clear(); 

	inputs_defined = false; coefficients_defined = false;
}

void DBR::set_params(double &min_feature, double &wl_Bragg, std::vector<double> &wl_vals, std::vector<double> &neff_w_low, std::vector<double> &neff_w_mid, std::vector<double> &neff_w_high, std::vector<double> &ngroup_w_mid)
{
	// Assign the input parameters to the DBR class
	// R. Sheehan 7 - 5 - 2019

	try {
		int n_wl = static_cast<int>(wl_vals.size()); 
		bool c1 = n_wl > 5 ? true : false;
		bool c2 = static_cast<int>(neff_w_low.size()) == n_wl ? true : false; 
		bool c3 = static_cast<int>(neff_w_mid.size()) == n_wl ? true : false; 
		bool c4 = static_cast<int>(neff_w_high.size()) == n_wl ? true : false; 
		bool c5 = static_cast<int>(ngroup_w_mid.size()) == n_wl ? true : false;
		bool c6 = wl_Bragg > wl_vals[0] && wl_Bragg < wl_vals[n_wl-1] ? true : false;
		bool c7 = min_feature > 0.0 ? true : false; 
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6 && c7; 

		if (c10) {

			n_wl_vals = n_wl; 

			// Assign the values of the inputs to the members
			Bragg_WL = wl_Bragg; 

			min_feat_size = min_feature; // minimum feature size possible during lithography

			// wavelength spectrum over which DBR can be calculated
			wavelengths = wl_vals; 

			// dispersion data for the different width waveguides computed over the wavelength spectrum
			neff_l = neff_w_low; 
			neff_m = neff_w_mid; 
			neff_h = neff_w_high; 

			// group index data for the mid width waveguide computed over the wavelength spectrum
			ngrp_m = ngroup_w_mid; 

			inputs_defined = true; 
		}
		else {
			std::string reason; 
			reason = "Error: void DBR::set_params(double &wl_Bragg, std::vector<double> &wl_vals, std::vector<double> &neff_w_low, std::vector<double> &neff_w_mid, std::vector<double> &neff_w_high)\n"; 
			if (!c1) reason += "n_wl: " + template_funcs::toString(n_wl) + " is not correct\n"; 
			if (!c2) reason += "neff_w_low.size(): " + template_funcs::toString(neff_w_low.size()) + " is not correct\n";
			if (!c3) reason += "neff_w_mid.size(): " + template_funcs::toString(neff_w_mid.size()) + " is not correct\n";
			if (!c4) reason += "neff_w_high.size(): " + template_funcs::toString(neff_w_high.size()) + " is not correct\n";
			if (!c5) reason += "ngroup_w_mid.size(): " + template_funcs::toString(ngroup_w_mid.size()) + " is not correct\n";
			if (!c6) reason += "wl_Bragg: " + template_funcs::toString(wl_Bragg,2) + " is outside range of input wavelength data\n";
			
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void DBR::coefficients()
{
	// compute the data for delta_n, delta_beta, kappa, gamma
	// this version will be easier to understand when I look back on it in a year's time
	
	// I could write a version that computes all coefficients in a single loop, this would be more efficient
	// but harder to understand when I look back on it later
	// in any case the computational overhead of multiple loops is bound to be negligible
	
	// R. Sheehan 8 - 5 - 2019

	try {
		if (inputs_defined) {

			index_contrast();

			beta_contrast();

			coupling_coeffs(); 

			gamma_coeffs(); 
			
			coefficients_defined = true;
		}
		else {
			std::string reason; 
			reason = "Error: void DBR::coefficients()\n";
			reason += "Inputs not correctly defined\n"; 

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void DBR::index_contrast()
{
	// compute the data for delta_n
	// R. Sheehan 5 - 4 - 2018

	try {
		if (inputs_defined) {

			double dn;
			for (int i = 0; i<n_wl_vals; i++) {
				dn = std::max(neff_h[i], neff_l[i]) - std::min(neff_h[i], neff_l[i]);
				dn /= 2.0*neff_m[i];
				delta_n.push_back(dn);
			}

		}
		else {
			std::string reason;
			reason = "Error: void DBR::compute_index_contrast()\n";
			reason += "Inputs not correctly defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void DBR::beta_contrast()
{
	// Compute delta beta for given Bragg WL
	// R. Sheehan 5 - 4 - 2018

	try {
		if (inputs_defined) {
			double delta_Bragg_neff, dB; 

			interpolation::polint(wavelengths, neff_m, Bragg_WL, Bragg_neff, delta_Bragg_neff); // compute neff at the Bragg WL

			Bragg_beta = beta(Bragg_WL, Bragg_neff); // beta_{0}

			period(); // compute the Bragg period and the grating order

			for (int i = 0; i < n_wl_vals; i++) {
				dB = beta(wavelengths[i], neff_m[i]) - Bragg_beta; 
				delta_beta.push_back(dB); 
			}
		}
		else {
			std::string reason;
			reason = "Error: void DBR::beta_contrast()\n";
			reason += "Inputs not correctly defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void DBR::coupling_coeffs()
{
	// Compute kappa from index contrast and grating type
	// delta_n contains the index constrast for waveguide
	// R. Sheehan 5 - 4 - 2018

	try {
		if (inputs_defined && static_cast<int>(delta_n.size()) == n_wl_vals ) {
			
			double kval; 
			for (int i = 0; i < n_wl_vals; i++) {
				
				kval = (2.0*delta_n[i]) / Bragg_WL; // coupling coefficient for straight sidewall corrugation
				
				kappa.push_back(kval); 
			}

		}
		else {
			std::string reason;
			reason = "Error: void DBR::coupling_coeffs()\n";
			reason += "Inputs not correctly defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void DBR::gamma_coeffs()
{
	// Compute gamma from kappa and delta beta
	// kappa is the coupling coefficient for corrugation type DBR_type 
	// Delta_beta is the propagation constant variation
	// in waveguide with Bragg wavelength WL_bragg
	// in waveguide with width WG_width and grating widths WG_h, WG_l at polarision pol
	// R. Sheehan 5 - 4 - 2018

	try {
		if (inputs_defined && static_cast<int>(kappa.size()) == n_wl_vals && static_cast<int>(delta_beta.size()) == n_wl_vals) {

			std::complex<double> gval;

			for (int i = 0; i < n_wl_vals; i++) {
				gval = template_funcs::DSQR(kappa[i]) - template_funcs::DSQR(delta_beta[i]);
				gval = std::sqrt(gval);

				gamma.push_back(gval);
			}
		}
		else {
			std::string reason;
			reason = "Error: void DBR::gamma_coeffs()\n";
			reason += "Inputs not correctly defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void DBR::r_spectrum(double &L_DBR)
{
	// compute the DBR reflectivity, actual value depends on number of grating periods used
	// R. Sheehan 8 - 5 - 2019

	try {
		if (coefficients_defined && L_DBR > 0.0 && Bragg_Period > min_feat_size) {

			n_periods = static_cast<int>( std::floor( L_DBR / Bragg_Period ) ); 

			L_DBR = n_periods * Bragg_Period; // adjust L_DBR to match the Bragg period and the no. periods

			std::complex<double> r;
			double absr;

			for (int i = 0; i < n_wl_vals; i++) {
				r = reflectivity(L_DBR, kappa[i], wavelengths[i], delta_beta[i], gamma[i]); 
				absr = std::abs(r); 
				r_spctrm.push_back( template_funcs::DSQR(absr) ); // Reflectance is given by |r|^{2}
			}
		}
		else {
			std::string reason;
			reason = "Error: void DBR::r_spectrum(double &L_DBR)\n";
			reason += "Inputs not correctly defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void DBR::Lambda_Rpeak_BW(double &L_DBR)
{
	// compute the DBR peak reflectivity and bandwidth
	// R. Sheehan 8 - 5 - 2019

	try {
		if (coefficients_defined && L_DBR > Bragg_Period && Bragg_Period > min_feat_size) {

			n_periods = static_cast<int>(std::floor(L_DBR / Bragg_Period));

			L_DBR = n_periods * Bragg_Period; // adjust L_DBR to match the Bragg period and the no. periods

			double kval, delta, ngval; 

			interpolation::polint(wavelengths, kappa, Bragg_WL, kval, delta); // compute kappa at WL_bragg from data

			interpolation::polint(wavelengths, ngrp_m, Bragg_WL, ngval, delta); // compute n_group at WL_bragg from data

			Rpeak = template_funcs::DSQR( std::tanh(kval*L_DBR) );

			BW = ( template_funcs::DSQR(Bragg_WL) / (PI*ngval) ) * sqrt( template_funcs::DSQR(kval) + template_funcs::DSQR(PI / L_DBR) );
		}
		else {
			std::string reason;
			reason = "Error: void DBR::Lambda_Rpeak_BW(double &L_DBR, double &Rpeak, double &BW)\n";
			reason += "Inputs not correctly defined\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void DBR::period()
{
	// compute the period of the grating structure from the computed parameters
	// also determine the order of the grating possible given the min_feat_size
	// R. Sheehan 8 - 5 - 2019

	try {
		if (inputs_defined && Bragg_neff > 0.0) {

			Bragg_Period = Bragg_WL / (2.0*Bragg_neff); // \Lambda
			Bragg_order = 1; 
			
			while (Bragg_Period < min_feat_size) {
				Bragg_Period *= 2.0; // Increase Bragg period so it's larger than min. feature size
				Bragg_order += 1; // increase order of grating
			}
		}
		else {
			std::string reason;
			reason = "Error: void DBR::period()\n";
			reason += "Bragg_neff: " + template_funcs::toString(Bragg_neff, 2) + " <= 0.0\n";
			throw std::runtime_error(reason);
		}
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what();
	}
}

double DBR::beta(double &lambda, double &neff)
{
	// compute propagation constant at Bragg wavelength in units of nm^{-1}
	// R. Sheehan 5 - 4 - 2018

	try {
		if (lambda > 0.0) {
			return ( ( Two_PI * neff ) / lambda ); 
		}
		else {
			std::string reason; 
			reason = "Error: double DBR::beta(double &lambda, double &neff)\n"; 
			reason += "lambda: " + template_funcs::toString(lambda, 2) + " <= 0.0\n"; 
			throw std::runtime_error(reason); 
		}
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what(); 
		return 0.0;
	}
}

std::complex<double> DBR::reflectivity(double &L_DBR, double &kval, double &wavel, double &delta_b, std::complex<double> &gval)
{
	// Compute complex reflectivity based on input data values
	// R. Sheehan 5 - 4 - 2018

	try{
		if (L_DBR > Bragg_Period) {
			std::complex<double> eye(0.0, 1.0);
			std::complex<double> sval, cval, gL, numer, denom;

			gL = gval * L_DBR;
			sval = sinh(gL); cval = cosh(gL);

			numer = -1.0*eye*kval*sval;
			denom = (gval*cval) + (eye*delta_b*sval);

			return numer / denom;
		}
		else {
			std::string reason; 
			reason = "Error: std::complex<double> DBR::reflectivity(double &L_DBR, double &kval, double &wavel, double &delta_b, std::complex<double> &gval)\n"; 
			reason += "L_DBR: " + template_funcs::toString(L_DBR, 2) + " <= 0.0\n"; 
			throw std::runtime_error(reason); 
		}
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what();
		return std::complex<double>(0.0, 0.0);
	}
}

void DBR::report()
{
	// Report on the results of a calculation
	// R. Sheehan 8 - 5 - 2019

	std::cout << "DBR Parameters\n\n"; 
	std::cout << "Grating WL: " << Bragg_WL << " nm\n"; 
	std::cout << "Grating Order: " << Bragg_order << "\n"; 
	std::cout << "Grating Period: " << Bragg_Period << " nm\n"; 
	std::cout << "Grating Count: " << n_periods << "\n";
	std::cout << "Grating Length: " << (n_periods * Bragg_Period)/1000.0 << " um\n"; 
	std::cout << "Grating BW: " << BW << " nm\n"; 
	std::cout << "Grating R: " << Rpeak << "\n\n"; 
}

void DBR::save_spctrm_to_file(std::string &filename)
{
	// save the computed reflectance spectrum to a file
	// R. Sheehan 8 - 5 - 2019

	try {
		bool c1 = wavelengths.size() == r_spctrm.size() ? true : false;
		bool c2 = filename != empty_str ? true : false;
		bool c3 = useful_funcs::valid_filename_length(filename);
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

			if (write.is_open()) {
				for (int i = 0; i < n_wl_vals; i++) {
					write << std::setprecision(10) << wavelengths[i] << " , " << r_spctrm[i] << "\n"; 
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
			if (!c2) reason += "filename not defined\n";
			if (!c3) reason += "filename: " + filename + " too long\n";
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