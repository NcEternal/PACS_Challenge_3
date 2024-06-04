#include "Parameter_Communication.hpp"
#include <mpi.h>

namespace laplace {

/* Function to let rank 0 communicate the parameters of the problem to the other ranks */
void communicate(parameters& p, int rank, int size) {
	
	/* No process to communicate to */
	if (size < 2)
		return;

	/* Broadcast domain bounds */
	MPI_Bcast(&p.x_lb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&p.x_ub, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&p.y_lb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&p.y_ub, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/* Broadcast mesh refinements */
	MPI_Bcast(&p.nx, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&p.ny, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	/* Size of the strings to communicate*/
	unsigned function_size;

	/* Communicate forcing term*/
	if (rank == 0)
		function_size = p.force.length();

	MPI_Bcast(&function_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	if (rank != 0)
		p.force.resize(function_size);

	MPI_Bcast(p.force.data(), function_size, MPI_CHAR, 0, MPI_COMM_WORLD);

	/* Communicate boundary terms */
	for (unsigned i = 0; i < 4; ++i) {
		if (rank == 0)
			function_size = p.bcs[i].length();

		MPI_Bcast(&function_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

		if (rank != 0)
			p.bcs[i].resize(function_size);

		MPI_Bcast(p.bcs[i].data(), function_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	}

	/* Broadcast boundary types */
	MPI_Bcast(p.bc_types.data(), 4, MPI_CHAR, 0, MPI_COMM_WORLD);

	/* Broadcast maximum iterations and tolerance */
	MPI_Bcast(&p.max_it, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&p.tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/* Broadcast exact solution */
	if (rank == 0)
		function_size = p.exact_sol.length();

	MPI_Bcast(&function_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	if (rank != 0)
		p.exact_sol.resize(function_size);

	MPI_Bcast(p.exact_sol.data(), function_size, MPI_CHAR, 0, MPI_COMM_WORLD);

}

/* Overload that gets rank and size on its own */
void communicate(parameters& p) {

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	communicate(p, rank, size);

}


}