#include "Solver.hpp"
#include "Parsers.hpp"
#include "Input_Output.hpp"
#include "Parameter_Communication.hpp"

#include <iostream>	// input-output

#include <mpi.h>	// parallelism

#include <numbers>	// pi
#include <cmath>	// sin

/* Main way to use the solver */

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/* Functions in the standard case, avoid using the parser */

	auto pi = std::numbers::pi;
	auto standard_f = [pi](double x, double y) { return 8 * pi * pi * std::sin(2 * pi * x) * std::sin(2 * pi * y); };
	std::function<double(double, double)> u_ex = [pi](double x, double y) { return std::sin(2 * pi * x) * std::sin(2 * pi * y); };
	auto standard_bc = [](double x) { return 0.; };

	std::function<double(double, double)> force = standard_f;
	std::vector<std::function<double(double)>> bcs({ standard_bc, standard_bc, standard_bc, standard_bc });



	/* Read Parameters */
	
	std::string filename{};

	if (rank == 0) {
		std::cout << "\n\n---- Laplace Equation Solver ----\n\n";
		std::cout << "Provide name of a parameter file in the './parameters/' folder (Leave empty for default problem): ";
		std::getline(std::cin, filename);
		std::cout << std::endl;
	}

	if (!filename.empty())
		filename = "parameters/" + filename;

	laplace::parameters default_;
	auto p = laplace::read(filename, rank);



	/* Communication */

	laplace::communicate(p);



	/* Check if functions were changed */

	bool function_change = false, solution_provided = false;

	if (p.force != default_.force) {
		function_change = true;
		force = parsing::Force(p.force);
	}
	for (unsigned i = 0; i < 4; ++i) {
		if (p.bcs[i] != default_.bcs[i]) {
			function_change = true;
			bcs[i] = parsing::Boundary_Condition(p.bcs[i], (i % 2 == 0) ? 'x' : 'y');
		}
	}
	if (p.bc_types != default_.bc_types)
		function_change = true;
	if (p.exact_sol != default_.exact_sol) {
		solution_provided = true;
		u_ex = parsing::Force(p.exact_sol);
	}



	/* Show Data */

	laplace::print_data(p, rank);



	/* Solve */

	const laplace::solver L(p.x_lb, p.x_ub, p.y_lb, p.y_ub, p.nx, p.ny, force, bcs[0], bcs[1], bcs[2], bcs[3], p.bc_types, rank, size);

	if (size > 1) 
		L.parallel_solve(p.max_it, p.tol);
	else 
		L.sequential_solve(p.max_it, p.tol);



	/* Results */

	if (rank == 0)
		std::cout << "Convergenge: " << L.has_converged();
	if (rank == 0 && (solution_provided == function_change))
		std::cout << "\nL2 Error: "	<< L.L2_error(u_ex) << std::endl;

	if (rank == 0) {
		std::string out_file;
		if (!filename.empty())
			out_file = "./output/out_" + filename.substr(11, filename.find('.', 11) - 11) + ".vtk";
		else
			out_file = "./output/out.vtk";
		std::cout << "\nSaving solution to \"" << out_file << "\"\n" << std::endl;
		L.to_VTK(out_file);
	}

	MPI_Finalize();

	return 0;

}
