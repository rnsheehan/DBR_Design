#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions for the slab base class

slab::slab(void)
{
	// Default Constructor
	M = 0;

	k_sqr = nc_sqr = ncl_sqr = ns_sqr = nr_sqr = nm_sqr = k_sqr_nc_sqr = k_sqr_ns_sqr = k_sqr_ncl_sqr = k_sqr_nr_sqr = 0.0;
	k_sqr_nm_sqr = aa = bb = d = dr = l = nc = ns = ncl = nr = nm = etacr = etacs = etarcl = 0.0;
	V = na = k = lower = upper = upper = w = efieldint = hfieldint = 0.0;
}

slab::~slab(void)
{
	// Deconstructor
	betaE.clear();

	betaH.clear(); 
}

int slab::nbeta(bool mode)
{
	// return the number of computed propagation constants for a given polarisation
	// it can be the case that the predicted value M is greater than the actual number of modes in a waveguide

	return static_cast<int>(mode ? betaE.size() : betaH.size() );
}

double slab::h(int i, bool t)
{
	//Wavenumber in Core

	try {

		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_tl_mode::h(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(k_sqr_nc_sqr - template_funcs::DSQR(betaE[i]));
				}
				else {
					return sqrt(k_sqr_nc_sqr - template_funcs::DSQR(betaH[i]));
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab::h(int i, bool t)\nNo modes have been computed\n");
			return 0; 
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::p(int i, bool t)
{
	//Wavenumber in Substrate

	try {
		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab::p(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(template_funcs::DSQR(betaE[i]) - k_sqr_ns_sqr);
				}
				else {
					return sqrt(template_funcs::DSQR(betaH[i]) - k_sqr_ns_sqr);
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab::p(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::q(int i, bool t)
{
	//Wavenumber in Cladding

	try {

		if (nbeta(t) > 0) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab::q(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt(template_funcs::DSQR(betaE[i]) - k_sqr_ncl_sqr);
				}
				else {
					return sqrt(template_funcs::DSQR(betaH[i]) - k_sqr_ncl_sqr);
				}
			}

		}
		else {
			throw std::invalid_argument("Error: double slab::q(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::r(int i, bool t)
{
	//Wavenumber in Rib

	try {

		if (nbeta(t) > 0) {

			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab::r(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return sqrt( template_funcs::DSQR(betaE[i]) - k_sqr_nr_sqr);
				}
				else {
					return sqrt( template_funcs::DSQR(betaH[i]) - k_sqr_nr_sqr);
				}
			}

		}
		else {
			throw std::invalid_argument("Error: double slab::r(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab::_beta(int i, bool t)
{
	// return the i^{th} propagation constant for a given polarisation
	// by convention t == true => TE modes, t == false => TM modes

	try {
		if (nbeta(t) > 0) {
			if (i<0 || i > nbeta(t)) {
				throw std::range_error("Error: double slab_wg::_beta(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (t) {
					return betaE[i];
				}
				else {
					return betaH[i];
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab::_beta(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Definitions for the three layer slab derived class
slab_tl_neff::slab_tl_neff(void)
{
	// Default Constructor
	g = 0; 
}

slab_tl_neff::slab_tl_neff(double width, double lambda, double ncore, double nsub, double nclad)
{
	set_params(width, lambda, ncore, nsub, nclad);
}

void slab_tl_neff::set_params(double width, double lambda, double ncore, double nsub, double nclad)
{
	// assign values to the parameters for the three layer slab waveguide
	// throw exception if one of the inputs is not correct

	try {

		bool c1 = width > 0.0 ? true : false;
		bool c2 = lambda > 0.0 ? true : false;
		bool c3 = nclad >= 1.0 ? true : false;
		bool c4 = nsub >= 1.0 ? true : false;
		bool c5 = ncore > std::max(nsub, nclad) ? true : false;

		if (c1 && c2 && c3 && c4 && c5) {

			d = width;

			l = lambda;

			nc = ncore;
			nc_sqr = template_funcs::DSQR(nc);

			ns = std::max(nsub, nclad);
			ns_sqr = template_funcs::DSQR(ns);

			ncl = std::min(nsub, nclad);
			ncl_sqr = template_funcs::DSQR(ncl);

			aa = (nc_sqr / ns_sqr);

			bb = (nc_sqr / ncl_sqr);

			k = Two_PI / l;
			k_sqr = template_funcs::DSQR(k); // k_{0}^{2}

			k_sqr_nc_sqr = k_sqr * nc_sqr; // k_{0}^{2} n_{c}^{2}
			k_sqr_ns_sqr = k_sqr * ns_sqr; // k_{0}^{2} n_{s}^{2}
			k_sqr_ncl_sqr = k_sqr * ncl_sqr; // k_{0}^{2} n_{cl}^{2}

			double x = nc_sqr - ns_sqr;
			double y = ns_sqr - ncl_sqr;

			// WG asymmetry paramater
			g = y > 0.0 ? (y / x) : 0.0;

			na = sqrt(x); // numerical aperture

			V = (PI*d*na) / l; // V-parameter

			// predicted number of modes
			M = y > 0 ? static_cast<int>( std::max(1.0, ceil((2.0*V / PI) - atan(g) / PI)) ) : static_cast<int>(std::max(1.0, ceil((2.0*V / PI))));

			lower = k * ns; // lower bound of search space

			upper = k * nc; // upper bound of search space

			w = k * SPEED_OF_LIGHT;

			efieldint = 0.5*EPSILON*SPEED_OF_LIGHT;

			hfieldint = 0.5*MU*SPEED_OF_LIGHT;

			// Empty the std::vector each time a new instance of the class is called
			betaE.clear();
			betaH.clear();
		}
		else {
			std::string reason = "Error: void slab_tl_neff::set_params(double width,double lambda,double ncore,double nsub,double nclad)\n";
			if (!c1) reason += "WG width = " + template_funcs::toString(width, 3) + " is negative\n";
			if (!c2) reason += "Wavelength = " + template_funcs::toString(lambda, 3) + " is negative\n";
			if (!c3) reason += "Cladding Index = " + template_funcs::toString(nclad, 3) + " is less than one\n";
			if (!c4) reason += "Substrate Index = " + template_funcs::toString(nsub, 3) + " is less than one\n";
			if (!c5) reason += "Core Index = " + template_funcs::toString(ncore, 3) + " is less than cladding / substrate index\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl_neff::eigeneqn_3(double x, bool t, int mm)
{
	// Using this version of the dispersion equation means that you don't have to do a bracketing step
	// which you would have to do if you used eigeneqn

	try {
		
		if (k_sqr_nc_sqr > k_sqr_ns_sqr) {

			double h, p, q, tmp;

			double x_sqr = template_funcs::DSQR(x);

			tmp = k_sqr_nc_sqr - x_sqr;
			h = (tmp>0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ns_sqr;
			p = (tmp>0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ncl_sqr;
			q = (tmp>0 ? sqrt(tmp) : 0.0);

			if (t) {//TE modes
				return d * h - mm * PI - atan(q / h) - atan(p / h);
			}
			else {//TM modes
				return d * h - mm * PI - atan(bb*(q / h)) - atan(aa*(p / h));
			}
		}
		else {
			std::string reason = "Error: double slab_tl_neff::eigeneqn_3(double x,bool t,int mm)\n";
			reason += "Input parameters not correctly defined\n";
			reason += "k_{0}^{2} n_{c}^{2} = " + template_funcs::toString(k_sqr_nc_sqr) + ", k_{0}^{2} n_{s}^{2} = " + template_funcs::toString(k_sqr_ns_sqr) + "\n";
			return 0;
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_tl_neff::zbrent(double x1, double x2, double tol, bool t, int mm)
{
	//Find the root of the function between x1 and x2 using brent's method
	//The root is refined to +/- tol
	//Seemingly one of the best methods to use

	//This will be used to compute the roots of eigeneq_3
	//R. Sheehan 28 - 5 - 2010

	try {

		bool c1 = fabs(x1 - x2) > 1.0e-9 ? true : false; // cannot have x1 == x2
		bool c2 = tol > 1.0e-16 ? true : false;
		bool c3 = mm < M ? true : false;

		if (c1 && c2 && c3) {

			int iter;

			static const int ITMAX = 100;//Used in rtbis, rtflsp, rtsec, zriddr

			double a = std::min(x1, x2), b = std::max(x1, x2), c = std::max(x1, x2), d, e, min1, min2;
			double fc, p, q, r, s, tol1, xm;
			double fa = eigeneqn_3(a, t, mm), fb = eigeneqn_3(b, t, mm);

			if ((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) {
				//std::cerr << "Root must be bracketed in zbrent\n";
			}
			fc = fb;
			for (iter = 1; iter <= ITMAX; iter++) {
				if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
					c = a;
					fc = fa;
					e = d = b - a;
				}
				if (fabs(fc)<fabs(fb)) {
					a = b;
					b = c;
					c = a;
					fa = fb;
					fb = fc;
					fc = fa;
				}
				tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
				xm = 0.5*(c - b);
				if (fabs(xm) <= tol1 || fb == 0.0) return b;
				/*if(fabs(xm)<=tol1 || fb==0.0){
				std::cout<<"Brent's Method converged in "<<iter<<" iterations\n";
				return b;
				}*/
				if (fabs(e) >= tol1 && fabs(fa)>fabs(fb)) {
					s = fb / fa;
					if (a == c) {
						p = 2.0*xm*s;
						q = 1.0 - s;
					}
					else {
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
						q = (q - 1.0)*(r - 1.0)*(s - 1.0);
					}
					if (p>0.0) q = -q;
					p = fabs(p);
					min1 = 3.0*xm*q - fabs(tol1*q);
					min2 = fabs(e*q);
					if (2.0*p<std::min(min1, min2)) {
						e = d;
						d = p / q;
					}
					else {
						d = xm;
						e = d;
					}
				}
				else {
					d = xm;
					e = d;
				}
				a = b;
				fa = fb;
				if (fabs(d)>tol1) {
					b += d;
				}
				else {
					b += template_funcs::SIGN(tol1, xm);
				}
				fb = eigeneqn_3(b, t, mm);
			}
			std::cerr << "Maximum number of iterations exceeded in zbrent\n";
			return 0.0;
		}
		else {
			std::string reason = "Error: double slab_tl_neff::zbrent(double x1,double x2,double tol,bool t,int mm)\n";
			if (!c1) reason += "Cannot have x1 = x2\nx1 = " + template_funcs::toString(x1) + ", x2 = " + template_funcs::toString(x1) + "\n";
			if (!c2) reason += "Desired tolerance is less than smallest allowable EPS\n";
			if (!c3) reason += "mm >= M\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_tl_neff::neff_search(bool mode)
{
	//This version solves the standard slab dispersion equation with a superior root finder
	//sic. Brent's Method NRinC ch. 9
	//R. Sheehan 28 - 5 - 2010

	try {

		if (lower < upper) {

			int m;
			double b;

			std::vector<double> vec;

			for (m = 0; m < M; m++) {
				b = zbrent(lower, upper, EPS, mode, m); 
				
				if (b>lower && b<upper) {
					vec.push_back(b);
				}

			}

			if (mode) {
				betaE = vec;
			}
			else {
				betaH = vec;
			}

			vec.clear();
		}
		else {
			std::string reason = "Error: void slab_tl_neff::neff_search(bool mode)\n";
			reason += "Search range is not correctly defined\n";
			reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n";
			throw std::range_error(reason);
		}
	}
	catch (std::range_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_tl_neff::report(bool mode)
{
	// print the computed results to screen
	std::cout << "Output for the slab waveguide calculator\n";
	std::cout << "\n";
	std::cout << "Waveguide width = " << d << " microns\n";
	std::cout << "Wavelength = " << l << " microns\n";
	std::cout << "Wavenumber = " << k << " inverse microns\n";
	std::cout << "\n";
	std::cout << "Core Index = " << nc << "\n";
	std::cout << "Substrate Index = " << ns << "\n";
	std::cout << "Cladding Index = " << ncl << "\n";
	std::cout << "\n";
	std::cout << "Numerical Aperture = " << na << "\n";
	std::cout << "V-Parameter = " << V << "\n";
	std::cout << "Number of Modes = " << M << "\n";
	std::cout << "\n";

	std::string pol = (mode ? "TE" : "TM");

	std::cout << pol << " Modes\n";
	std::cout << "There are " << nbeta(mode) << " calculated " + pol + " modes\n";

	for (int i = 0; i<nbeta(mode); i++) {
		std::cout << pol + "_{" << i + 1 << "} = " << std::setprecision(10) << _beta(i, mode) << " , n_eff_{" << i + 1 << "} = " << std::setprecision(10) << (_beta(i, mode) / k) << "\n";
	}
	std::cout << "\n";
}

int slab_tl_neff::get_nmodes(bool mode)
{
	// return the number of computed propagation constants for a given polarisation
	// it can be the case that the predicted value M is greater than the actual number of modes in a waveguide

	return nbeta(mode); 
}

double slab_tl_neff::get_neff(int i, bool mode)
{
	// return the i^{th} effective index for a given polarisation
	// by convention mode == true => TE modes, mode == false => TM modes

	try {
		if (nbeta(mode) > 0) {
			if (i<0 || i > nbeta(mode)) {
				throw std::range_error("Error: double slab_tl_neff::get_neff(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (mode) {
					return betaE[i] / k;
				}
				else {
					return betaH[i] / k;
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_tl_neff::get_neff(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// slab_fl_mode_B is used to compute the effective indices in the four layer slab case B

slab_fl_neff_B::slab_fl_neff_B()
{
	// Default Constructor
}

slab_fl_neff_B::slab_fl_neff_B(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib)
{
	// Constructor

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)

	set_params(width, rib_width, lambda, ncore, nsub, nclad, nrib);
}

void slab_fl_neff_B::set_params(double width, double rib_width, double lambda, double ncore, double nsub, double nclad, double nrib)
{
	// set parameters for four layer slab

	// Case A => Field Oscillating in Core and Ridge
	// For there to be a solution one has to have ns <= ncl < nr < nc

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)

	try {
		bool c1 = width > 0.0 ? true : false;
		bool c1a = rib_width > 0.0 ? true : false;
		bool c2 = lambda > 0.0 ? true : false;
		bool c3 = nclad >= 1.0 ? true : false;
		bool c4 = nsub >= 1.0 ? true : false;
		bool c6 = nrib >= 1.0 ? true : false;
		bool c5 = ncore > std::max(nsub, nclad) ? true : false;
		bool c7 = ncore > std::max(nsub, nrib) ? true : false;

		if (c1 && c1a && c2 && c3 && c4 && c5 && c6 && c7) {

			d = width;

			dr = rib_width;

			l = lambda;

			nc = ncore;
			nc_sqr = template_funcs::DSQR(nc);

			nr = nrib;
			nr_sqr = template_funcs::DSQR(nr);

			ns = std::max(nsub, nclad);
			ns_sqr = template_funcs::DSQR(ns);

			ncl = std::min(nsub, nclad);
			ncl_sqr = template_funcs::DSQR(ncl);

			nm = std::max(ns, nr);
			nm_sqr = template_funcs::DSQR(nm);

			k = Two_PI / l;
			k_sqr = template_funcs::DSQR(k); // k_{0}^{2}

			k_sqr_nc_sqr = k_sqr * nc_sqr; // k_{0}^{2} n_{c}^{2}
			k_sqr_ns_sqr = k_sqr * ns_sqr; // k_{0}^{2} n_{s}^{2}
			k_sqr_ncl_sqr = k_sqr * ncl_sqr; // k_{0}^{2} n_{cl}^{2}
			k_sqr_nr_sqr = k_sqr * nr_sqr; // k_{0}^{2} n_{r}^{2}
			k_sqr_nm_sqr = k_sqr * nm_sqr; // k_{0}^{2} n_{m}^{2}

			etacs = nc_sqr / ns_sqr;
			etacr = nc_sqr / nr_sqr;
			etarcl = nr_sqr / ncl_sqr;

			// Only difference between A, B cases is search space for neff and eigenequation
			// Case A: lower = k ncl, upper = k nr, NA^{2} = nr^{2} - ncl^{2}
			// Case B: lower = k nm, upper - k nc, NA^{2} = nc^{2} - nm^{2}

			double x = nc_sqr - nm_sqr;
			double y = ns_sqr - ncl_sqr;

			na = sqrt(x); // numerical aperture

			V = (PI*d*na) / l; // V-parameter

			// predicted number of modes
			M = static_cast<int>(std::max(1.0, ceil((2.0*V / PI))));

			lower = k * nm; // lower bound of search space

			upper = k * nc; // upper bound of search space

			w = k * SPEED_OF_LIGHT;

			efieldint = 0.5*EPSILON*SPEED_OF_LIGHT;

			hfieldint = 0.5*MU*SPEED_OF_LIGHT;

			// Empty the std::vector each time a new instance of the class is called
			betaE.clear();
			betaH.clear();
		}
		else {
			std::string reason = "Error: void slab_fl_neff_B::set_params(double width,double lambda,double ncore,double nsub,double nclad)\n";
			if (!c1) reason += "WG width = " + template_funcs::toString(width, 3) + " is negative\n";
			if (!c2) reason += "Wavelength = " + template_funcs::toString(lambda, 3) + " is negative\n";
			if (!c3) reason += "Cladding Index = " + template_funcs::toString(nclad, 3) + " is less than one\n";
			if (!c4) reason += "Substrate Index = " + template_funcs::toString(nsub, 3) + " is less than one\n";
			if (!c6) reason += "Rib Index = " + template_funcs::toString(nrib, 3) + " is less than one\n";
			if (!c5) reason += "Core Index = " + template_funcs::toString(ncore, 3) + " is less than cladding / substrate index\n";
			if (!c7) reason += "Core Index = " + template_funcs::toString(ncore, 3) + " is less than rib / substrate index\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_fl_neff_B::eigeneqn_b(double x, int mm, bool t)
{
	//Dispersion equation corresponding to case b from Adams
	//This means that the field oscillates in the core only	
	//This is an alternative form that produces the correct solutions

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)

	try {

		if (k_sqr_nc_sqr > k_sqr_nm_sqr) {
			double h, p, q, r, tmp, x_sqr, v, v1, v2den, v2;

			x_sqr = template_funcs::DSQR(x);

			tmp = k_sqr_nc_sqr - x_sqr;
			h = (tmp > 0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ns_sqr;
			p = (tmp > 0 ? sqrt(tmp) : 0.0);

			tmp = x_sqr - k_sqr_ncl_sqr;
			q = (tmp > 0 ? sqrt(tmp) : 0.0);
			//if (!t) q *= etarcl; // multiply q by etarcl in the case of TM polarisation

			tmp = x_sqr - k_sqr_nr_sqr;
			r = (tmp > 0 ? sqrt(tmp) : 0.0);

			v = ((r - q) / (r + q)); // this includes the change to q depending on polarisation

			v1 = exp(-2.0 * r *dr); // In the limit of large dr v1 -> 0 and v2 -> 1 => Case B FLS -> TLS

			v2den = (1 + v * v1);

			v2 = (v2den > 0.0 ? (1 - v * v1) / v2den : 10.0);

			if (t) {//TE modes
				return (d * h) - atan((p / h)) - atan((r / h)*v2) - mm * PI;
			}
			else {//TM modes
				return (d * h) - atan(etacs*(p / h)) - atan(etacr * (r / h) * v2) - mm * PI;
			}
		}
		else {
			std::string reason = "Error: double slab_fl_neff_B::eigeneqn_b(double x, int mm, bool t)\n";
			reason += "Input parameters not correctly defined\n";
			reason += "k_{0}^{2} n_{c}^{2} = " + template_funcs::toString(k_sqr_nc_sqr) + ", k_{0}^{2} n_{m}^{2} = " + template_funcs::toString(k_sqr_nm_sqr) + "\n";
			return 0;
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double slab_fl_neff_B::zbrent(double x1, double x2, double tol, bool t, int mm)
{
	//Find the root of the function between x1 and x2 using brent's method
	//The root is refined to +/- tol
	//Seemingly one of the best methods to use

	//This will be used to compute the roots of eigeneq_b
	//R. Sheehan 28 - 5 - 2010

	try {

		bool c1 = fabs(x1 - x2) > 1.0e-9 ? true : false; // cannot have x1 == x2
		bool c2 = tol > 1.0e-16 ? true : false;
		bool c3 = mm < M ? true : false;

		if (c1 && c2 && c3) {

			int iter;

			static const int ITMAX = 100;//Used in rtbis, rtflsp, rtsec, zriddr

			double a = std::min(x1, x2), b = std::max(x1, x2), c = std::max(x1, x2), d, e, min1, min2;
			double fc, p, q, r, s, tol1, xm;
			double fa = eigeneqn_b(a, t, mm), fb = eigeneqn_b(b, t, mm);

			if ((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) {
				std::cerr << "Root must be bracketed in zbrent\n";
			}
			fc = fb;
			for (iter = 1; iter <= ITMAX; iter++) {
				if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
					c = a;
					fc = fa;
					e = d = b - a;
				}
				if (fabs(fc)<fabs(fb)) {
					a = b;
					b = c;
					c = a;
					fa = fb;
					fb = fc;
					fc = fa;
				}
				tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
				xm = 0.5*(c - b);
				if (fabs(xm) <= tol1 || fb == 0.0) return b;
				/*if(fabs(xm)<=tol1 || fb==0.0){
				std::cout<<"Brent's Method converged in "<<iter<<" iterations\n";
				return b;
				}*/
				if (fabs(e) >= tol1 && fabs(fa)>fabs(fb)) {
					s = fb / fa;
					if (a == c) {
						p = 2.0*xm*s;
						q = 1.0 - s;
					}
					else {
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
						q = (q - 1.0)*(r - 1.0)*(s - 1.0);
					}
					if (p>0.0) q = -q;
					p = fabs(p);
					min1 = 3.0*xm*q - fabs(tol1*q);
					min2 = fabs(e*q);
					if (2.0*p<std::min(min1, min2)) {
						e = d;
						d = p / q;
					}
					else {
						d = xm;
						e = d;
					}
				}
				else {
					d = xm;
					e = d;
				}
				a = b;
				fa = fb;
				if (fabs(d)>tol1) {
					b += d;
				}
				else {
					b += template_funcs::SIGN(tol1, xm);
				}
				fb = eigeneqn_b(b, t, mm);
			}
			std::cerr << "Maximum number of iterations exceeded in zbrent\n";
			return 0.0;
		}
		else {
			std::string reason = "Error: double slab_fl_neff_B::zbrent(double x1,double x2,double tol,bool t,int mm)\n";
			if (!c1) reason += "Cannot have x1 = x2\nx1 = " + template_funcs::toString(x1) + ", x2 = " + template_funcs::toString(x1) + "\n";
			if (!c2) reason += "Desired tolerance is less than smallest allowable EPS\n";
			if (!c3) reason += "mm >= M\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_fl_neff_B::neff_search(bool mode)
{
	// Compute the waveguide mode effective indices based on given polarisation and wg_type

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)

	// Case B: lower = k nm, upper - k nc, NA^{2} = nc^{2} - nm^{2}

	try {

		if (lower < upper) {

			int m;
			double b;

			std::vector<double> vec;

			for (m = 0; m < M; m++) {
				b = zbrent(lower, upper, EPS, mode, m);

				if (b>lower && b<upper) {
					vec.push_back(b);
				}

			}

			if (mode) {
				betaE = vec;
			}
			else {
				betaH = vec;
			}

			vec.clear();

		}
		else {
			std::string reason = "Error: void slab_fl_neff_B::neff_search(bool mode, bool wg_type)\n";
			reason += "Search range is not correctly defined\n";
			reason += "lower = " + template_funcs::toString(lower, 4) + ", upper = " + template_funcs::toString(upper, 4) + "\n";
			throw std::range_error(reason);
		}
	}
	catch (std::range_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void slab_fl_neff_B::report(bool mode)
{
	// print the computed results to screen
	std::cout << "Output for the slab waveguide calculator\n";
	std::cout << "\n";
	std::cout << "Waveguide width = " << d << " microns\n";
	std::cout << "Rib width = " << dr << " microns\n";
	std::cout << "Wavelength = " << l << " microns\n";
	std::cout << "Wavenumber = " << k << " inverse microns\n";
	std::cout << "\n";
	std::cout << "Core Index = " << nc << "\n";
	std::cout << "Substrate Index = " << ns << "\n";
	std::cout << "Rib Index = " << nr << "\n";
	std::cout << "Cladding Index = " << ncl << "\n";
	std::cout << "\n";
	std::cout << "Numerical Aperture = " << na << "\n";
	std::cout << "V-Parameter = " << V << "\n";
	std::cout << "Number of Modes = " << M << "\n";
	std::cout << "\n";

	std::string pol = (mode ? "TE" : "TM");

	std::cout << pol << " Modes\n";
	std::cout << "There are " << nbeta(mode) << " calculated " + pol + " modes\n";

	for (int i = 0; i < nbeta(mode); i++) {
		std::cout << pol + "_{" << i + 1 << "} = " << std::setprecision(10) << _beta(i, mode) << " , n_eff_{" << i + 1 << "} = " << std::setprecision(10) << (_beta(i, mode) / k) << "\n";
	}
	std::cout << "\n";
}

int slab_fl_neff_B::get_nmodes(bool mode)
{
	// return the number of computed propagation constants for a given polarisation
	// it can be the case that the predicted value M is greater than the actual number of modes in a waveguide

	return nbeta(mode);
}

double slab_fl_neff_B::get_neff(int i, bool mode)
{
	// return the i^{th} effective index for a given polarisation
	// by convention mode == true => TE modes, mode == false => TM modes

	try {
		if (nbeta(mode) > 0) {
			if (i<0 || i > nbeta(mode)) {
				throw std::range_error("Error: double slab_fl_neff_B::get_neff(int i, bool t)\n Attempting to access arrays out of range\n");
				return 0;
			}
			else {
				if (mode) {
					return betaE[i] / k;
				}
				else {
					return betaH[i] / k;
				}
			}
		}
		else {
			throw std::invalid_argument("Error: double slab_fl_neff_B::get_neff(int i, bool t)\nNo modes have been computed\n");
			return 0;
		}
	}
	catch (std::range_error &e) {
		std::cerr << e.what();
		return 0;
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}