#pragma once
#ifndef HH_Parameters_HH
#define HH_Parameters_HH

#include <string>
#include <vector>

namespace laplace {

/* Struct to represent the data of a Laplace Problem */

struct parameters {

	double x_lb = 0.;											// Lower x bound of the domain
	double x_ub = 1.;											// Upper x bound of the domain
	double y_lb = 0.;											// Lower y bound of the domain
	double y_ub = 1.;											// Upper y bound of the domain
	unsigned nx = 64;											// Number of points on the x axis
	unsigned ny = 64;											// Number of points on the y axis
	std::string force{ "8*_pi^2*sin(2*_pi*x)*sin(2*_pi*y)" };	// Force
	std::vector<std::string> bcs = { "0", "0", "0", "0" };		// Boundary Conditions, given in clockwise order starting from the y = y_top boundary
	std::string bc_types{ "dddd" };								// Boundary types. Same order as bcs. 'd' = Dirichlet, 'n' = Neumann
	unsigned max_it = 1e5;										// Maximum number of iterations
	double tol = 1e-6;											// Tolerance for convergence check
	std::string exact_sol{ "sin(2*_pi*x)*sin(2*_pi*y)" };		// Exact Solution to calculate the error at the end of the solving process

};

}

#endif