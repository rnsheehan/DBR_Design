#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::disp_curve_wire()
{
	// Compute the dispersion curve for a wire waveguide object
	// R. Sheehan 7 - 3 - 2019

	bool pol = TM;

	int n_pts;
	double start, stop, W, H;

	sweep WL;
	wg_dims wire_dims;

	Air ri_air;
	SiN ri_si;
	SiO2 ri_sio2;

	n_pts = 50; start = 1.5; stop = 1.6;
	WL.set_vals(n_pts, start, stop);

	W = 2; H = 0.3;

	wire_dims.set_rect_wire(W, H);

	wire_dispersion disp_calc;

	std::string filename;

	filename = "Wire_WG_Dispersion_W_" + template_funcs::toString(W, 2) + "_H_" + template_funcs::toString(H, 2) + dottxt;

	disp_calc.compute_dispersion_data(pol, WL, wire_dims, &ri_si, &ri_sio2, &ri_air, filename);

	disp_calc.get_computed_data(); 

	disp_calc.clear(); 
}

void testing::disp_curve_rib()
{
	// Compute the dispersion curve for a rib waveguide object
	// R. Sheehan 8 - 3 - 2019

	bool pol = TM;

	int n_pts;
	double start, stop, W, E, T;

	sweep WL;
	wg_dims rib_dims;

	Air ri_air;
	Si ri_si;
	SiO2 ri_sio2;

	n_pts = 50; start = 1.4; stop = 1.6;
	WL.set_vals(n_pts, start, stop);

	W = 1; E = 0.5; T = 0.3;
	rib_dims.set_rib(W, E, T);

	rib_dispersion disp_calc;

	std::string filename = "Rib_WG_Dispersion.txt";

	disp_calc.compute_dispersion_data(pol, WL, rib_dims, &ri_si, &ri_sio2, &ri_air, filename);
}

void testing::DBR_wire()
{
	// Compute the DBR properties for a wire based DBR
	// R. Sheehan 8 - 5 - 2019

	bool pol = TE;

	int n_pts;
	double start, stop, W, H;

	sweep WL;
	wg_dims wire_dims;

	Air ri_air;
	SiN ri_si;
	SiO2 ri_sio2;

	wire_dispersion disp_calc;

	std::vector<std::vector<double>> the_data_l;
	std::vector<std::vector<double>> the_data_m;
	std::vector<std::vector<double>> the_data_h;

	std::string filename;

	// for dispersion calculation it is assumed that wavelength is in units of um
	n_pts = 300; start = 1.4; stop = 1.68;
	WL.set_vals(n_pts, start, stop);

	// First Dispersion Calculation
	W = 1; H = 0.3;

	wire_dims.set_rect_wire(W, H);	

	filename = "Wire_WG_Dispersion_W_" + template_funcs::toString(W, 2) + "_H_" + template_funcs::toString(H, 2) + dottxt;

	disp_calc.compute_dispersion_data(pol, WL, wire_dims, &ri_si, &ri_sio2, &ri_air, filename);
	
	the_data_m = disp_calc.get_computed_data();

	disp_calc.clear(); 

	// Second Dispersion Calculation
	W = 0.5; H = 0.3;

	wire_dims.set_rect_wire(W, H);

	filename = "Wire_WG_Dispersion_W_" + template_funcs::toString(W, 2) + "_H_" + template_funcs::toString(H, 2) + dottxt;

	disp_calc.compute_dispersion_data(pol, WL, wire_dims, &ri_si, &ri_sio2, &ri_air, filename);

	the_data_l = disp_calc.get_computed_data();

	disp_calc.clear();

	// Third Dispersion Calculation
	W = 1.5; H = 0.3;

	wire_dims.set_rect_wire(W, H);

	filename = "Wire_WG_Dispersion_W_" + template_funcs::toString(W, 2) + "_H_" + template_funcs::toString(H, 2) + dottxt;

	disp_calc.compute_dispersion_data(pol, WL, wire_dims, &ri_si, &ri_sio2, &ri_air, filename);

	the_data_h = disp_calc.get_computed_data();

	disp_calc.clear();

	// DBR Reflectivity Calculation
	DBR grating; 

	double min_x, wl_brg, L_dbr;

	min_x = 100; wl_brg = 1532; 
	grating.set_params(min_x, wl_brg, the_data_m[0], the_data_l[1], the_data_m[1], the_data_h[1], the_data_m[2]); 

	grating.coefficients(); 

	L_dbr = 18*1000; // value in um converted to nm

	grating.r_spectrum(L_dbr); 

	grating.Lambda_Rpeak_BW(L_dbr); 

	grating.report();

	filename = "Wire_WG_DBR_Refl" + dottxt;

	grating.save_spctrm_to_file(filename); 
}

void testing::DBR_rib()
{
	// Compute the DBR properties for a wire based DBR
	// R. Sheehan 8 - 5 - 2019

	bool pol = TM;

	int n_pts;
	double start, stop, W, E, T;

	sweep WL;
	wg_dims rib_dims;

	Air ri_air;
	Si ri_si;
	SiO2 ri_sio2;

	rib_dispersion disp_calc;

	std::vector<std::vector<double>> the_data_l;
	std::vector<std::vector<double>> the_data_m;
	std::vector<std::vector<double>> the_data_h;

	std::string filename;

	n_pts = 300; start = 1.4; stop = 1.68;
	WL.set_vals(n_pts, start, stop);

	// First Dispersion Calculation
	W = 0.5; E = 0.5; T = 0.3;
	rib_dims.set_rib(W, E, T);

	filename = "Rib_WG_Dispersion.txt";

	disp_calc.compute_dispersion_data(pol, WL, rib_dims, &ri_si, &ri_sio2, &ri_air, filename);

	the_data_l = disp_calc.get_computed_data();

	disp_calc.clear();

	// Second Dispersion Calculation
	W = 1.0; E = 0.5; T = 0.3;
	rib_dims.set_rib(W, E, T);

	filename = "Rib_WG_Dispersion.txt";

	disp_calc.compute_dispersion_data(pol, WL, rib_dims, &ri_si, &ri_sio2, &ri_air, filename);

	the_data_m = disp_calc.get_computed_data();

	disp_calc.clear();

	// Third Dispersion Calculation
	W = 2.0; E = 0.5; T = 0.3;
	rib_dims.set_rib(W, E, T);

	filename = "Rib_WG_Dispersion.txt";

	disp_calc.compute_dispersion_data(pol, WL, rib_dims, &ri_si, &ri_sio2, &ri_air, filename);

	the_data_h = disp_calc.get_computed_data();

	disp_calc.clear();

	// DBR Reflectivity Calculation
	DBR grating;

	double min_x, wl_brg, L_dbr;

	min_x = 100; wl_brg = 1532;
	grating.set_params(min_x, wl_brg, the_data_m[0], the_data_l[1], the_data_m[1], the_data_h[1], the_data_m[2]);

	grating.coefficients();

	L_dbr = 36 * 1000; // value in um converted to nm

	grating.r_spectrum(L_dbr);

	grating.Lambda_Rpeak_BW(L_dbr);

	grating.report();

	filename = "Rib_WG_DBR_Refl" + dottxt;

	grating.save_spctrm_to_file(filename);
}


