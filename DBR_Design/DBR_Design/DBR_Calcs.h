#ifndef DBR_CALCS_H
#define DBR_CALCS_H

// Definition of a class used to compute DBR parameters based on computed waveguide dispersion data
// Implementation is based on an older code I wrote, see files DBR_Analysis.h and DBR_Analysis.cpp from old EIM code
// Theory is based on the PhD Thesis of Xu Wang, "Si Photonic Waveguide Bragg Gratings", 2013
// There's no nice way to do the calculation from a single class
// I have to do the dispersion calc based on the wire_dispersion class or the rib_dispersion class
// It make sense to follow the same structure and do the DBR calculations using dynamic typing so DBR will be fully private
// The idea is to have a class for doing DBR based on wire waveguide and another for doing DBR based on rib waveguides
// This can be extended whenever another waveguide structure is required
// A single DBR class will do the DBR parameter calculations based on the computed dispersion data for either the rib or wire waveguide
// DBR will ony take in dispersion data as input and do all the calculations accordingly
// Derived classes in will perform dispersion calculations and pass the data to the dispersion class
// I could just do wire and rib DBR as methods in a namespace, there's a lot of needless passing of information around. 
// Six of one, half dozen of another I'd say. 
// R. Sheehan 7 - 5 - 2019

class DBR {
protected:
	DBR();

	~DBR(); 

	void set_params(double &wl_res, sweep &swp_obj, std::vector<double> &neff_w_low, std::vector<double> &neff_w_mid, std::vector<double> &neff_w_high);

	void compute_coefficients(); // compute the data for delta_, delta_beta, kappa, gamma

	void compute_r_spectrum(double &L_DBR); //  compute the DBR reflectivity, actual value depends on number of grating periods used

	void compute_Lambda_Rpeak_BW(double &L_DBR, double &Period, double &Rpeak, double &BW); // compute the DBR peak reflectivity and bandwidth


private:
	//double L; // length of DBR structure in units of nm
	double Bragg_WL; // resonance wavelength of DBR structure
	//double Period; // period of DBR structure in units of nm
	//double Rpeak; // peak reflectivity of DBR structure
	//double BW; // spectral bandwidth of DBR structure

	std::vector<double> wavelengths; // wavelengths over which effective index values are computed
	std::vector<double> neff_m; // effective indices of the main input waveguide
	std::vector<double> neff_l; // effective indices of the low width waveguide
	std::vector<double> neff_h; // effective indices of the high width waveguide
	
	std::vector<double> delta_n; // index contrast for the DBR structure
	std::vector<double> delta_beta; // propagation constant contrast for the DBR structure
	std::vector<double> kappa; // coupling coefficients for the DBR structure
	std::vector<double> r_spectrum; // reflectivity spectrum for the DBR structure

	std::vector<std::complex<double>> gamma; // coefficients used in the calculation of reflectivity spectrum

};

#endif
