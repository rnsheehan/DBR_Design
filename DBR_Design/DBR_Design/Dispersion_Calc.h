#ifndef DISPERSION_CALC_H
#define DISPERSION_CALC_H

// Declaration of the dispersion calculation object
// The aim is to compute the effective index and group index versus wavelength curves for a 
// given waveguide structure. Data is waveguide geometry and material dependent
// R. Sheehan 22 - 2 - 2019

// Update the EIM code to include virtual functions for reducing a wg and assigning parameters
// That way you can make EIM object part of the dispersion class as a pointer to an object derived from EIM
// It will make the coding much easier and obviate the need to have derived classes for the dispersion class
// R. Sheehan 22 - 2 - 2019
// Done R. Sheehan 28 - 2 - 2019

class dispersion {
public:
	dispersion(); 

	dispersion(sweep &swp_obj, material *Ncore, material *Nsub, material *Nclad);

	~dispersion(); 

	void set_params(sweep &swp_obj, material *Ncore, material *Nsub, material *Nclad);

	void compute_dispersion(bool polarisation, wg_dims &dim_obj, std::string &filename, bool loud = false);

	std::vector<std::vector<double>> get_computed_data(bool cnvrt_wl_nm = true); // return all computed dispersion data

	void clear(); 

private:
	void compute_group_index(bool loud = false); 

	void save_data_to_file(std::string &filename);

protected:
	bool params_defined; 

	// the sweep object defines the wavelength sweep space, all wavelength values are in units of um
	sweep wavelength; 

	material *core; // object for the core material
	material *substrate; // object for the substrate material
	material *cladding; // object for the cladding material

	EIM *neff_calc; // object for computing the waveguide effective index 

	std::vector<double> neff_vals; 
	std::vector<double> ng_vals; 
};

class wire_dispersion : public dispersion {
public:
	wire_dispersion(); 

	void compute_dispersion_data(bool polarisation, sweep &swp_obj, wg_dims &dim_obj, material *Ncore, material *Nsub, material *Nclad, std::string &filename, bool loud = false);
};

class rib_dispersion : public dispersion {
public:
	rib_dispersion();

	void compute_dispersion_data(bool polarisation, sweep &swp_obj, wg_dims &dim_obj, material *Ncore, material *Nsub, material *Nclad, std::string &filename, bool loud = false);

};

#endif
