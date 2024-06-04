#include "Solver.hpp"
#include <mpi.h>		// parallelism
#include <fstream>		// writing to a file
#include <iostream>		// cerr
#include <cmath>		// sqrt



namespace laplace{

/* Pre-calculates Force over the Mesh */
void solver::pre_calc_force() const {

	/* Initialize force matrix to be the same size as m_solution (Which is always at least default initialized)*/
	m_force_matrix = std::vector<double>(m_solution.size(), 0.);	
	/* Get mesh sizes */
	unsigned n_cols = m_mesh.nx(), n_rows = m_mesh.ny();

	/* Sequential initialization case */
	if (size < 2) {

		/* Calculate Force */
		for (unsigned j = 0; j < n_rows; ++j) 
			for (unsigned i = 0; i < n_cols; ++i) 
				m_force_matrix[j * n_cols + i] = m_force(m_mesh.xi(i), m_mesh.yj(j));

	}
	/* Parallel initialization case */
	else {

		/* Get Local number of rows */
		unsigned local_rows = (static_cast<unsigned>(rank) < (n_rows % size)) ? n_rows / size + 1 : n_rows / size;
		/* Offset due to the fact that ranks other than 0 have an extra row at the bottom (or top depending on how you look at the matrix, I think of it going up from 0 to n_rows) */
		unsigned not_rank0 = (rank != 0) ? 1 : 0;

		/* Get row_offset to calculate the correct value of the force */
		unsigned row_offset;
		if (rank == 0)
			row_offset = 0;
		else if (rank == size - 1)
			row_offset = n_rows - local_rows;
		else
			row_offset = (local_rows == (n_rows / size) + 1) ? (local_rows * rank) : ((local_rows + 1) * (n_rows % size) + local_rows * (rank - (n_rows % size)));

		/* Calculate Force */
		for (unsigned j = 0; j < local_rows; ++j)
			for (unsigned i = 0; i < n_cols; ++i)
				m_force_matrix[(j + not_rank0) * n_cols + i] = m_force(m_mesh.xi(i), m_mesh.yj(j + row_offset));

	}

}



/* Assigns Boundary Conditions */
void solver::assign_bc() const {

	/* Get mesh sizes */
	unsigned n_cols = m_mesh.nx(), n_rows = m_mesh.ny();

	/* Sequential initialization case */
	if (size < 2) {

		/* Size of solution = size of mesh */
		m_solution = std::vector<double>(n_rows * n_cols, 0.);

		/* All 4 boundaries must be checked */
		for (unsigned b : {0, 1, 2, 3})
			assign_bc_helper(b, n_rows, n_cols, 0, 0);

	}
	/* Parallel initialization case */
	else {

		/* Get local_rows */
		unsigned local_rows = (static_cast<unsigned>(rank) < (n_rows % size)) ? (n_rows / size) + 1 : (n_rows / size);

		/* Master Rank */
		if (rank == 0) {

			/* Size of solution = size of mesh so that it can then be collected more easily */
			m_solution = std::vector<double>(n_rows * n_cols, 0.);
		
			/* The top boundary doesn't need to be checked since it doesn't belong to rank 0 */
			for (unsigned b : {1, 2, 3})
				assign_bc_helper(b, local_rows, n_cols, 0, 0);

		}
		/* Last Rank */
		else if (rank == size - 1) {

			/* Row offset to get correct value of boundary condition */
			unsigned row_offset = n_rows - local_rows;
			/* m_solution is already initialized since all parallel constructors call parallel_init() first and there are no other ways to change its size
			   so this just resets it in case parallel_solve was called and the boundary conditions were changed */
			m_solution = std::vector<double>(m_solution.size(), 0.);

			/* Bottom boundary doesn't need to be checked since it doesn't belong to rank size - 1 */
			for (unsigned b : {0, 1, 3})
				assign_bc_helper(b, local_rows, n_cols, row_offset, 1);

		}
		/* General Ranks */
		else {

			/* Row offset to get correct value of boundary condition */
			unsigned row_offset = (local_rows == (n_rows / size + 1)) ? (local_rows * rank) : ((local_rows + 1) * (n_rows % size) + local_rows * (rank - (n_rows % size)));
			/* Same as last rank */
			m_solution = std::vector<double>(m_solution.size(), 0.);

			/* Only side boundaries need to be checked */
			for (unsigned b : {1, 3})
				assign_bc_helper(b, local_rows, n_cols, row_offset, 1);

		}
	}

}



/* Helper to Assign Boundary Conditions */
void solver::assign_bc_helper(const unsigned boundary, const unsigned n_rows, const unsigned n_cols, const unsigned row_offset, const unsigned not_rank0) const {

	/* The varying coordinate (default x) */
	unsigned varying = 0;
	/* The fixed coordinate (default = y) */
	unsigned fixed = (boundary > 1) ? 0 : n_rows - 1;
	/* The size of the varying coordinate */
	unsigned n = n_cols;
	/* A lambda to extract the correct coordinate value */
	std::function<double(unsigned)> coord = [this](unsigned i) { return m_mesh.xi(i); };

	/* If the boundary is odd, then invert the roles of x and y */
	if (boundary % 2 == 1) {
		n = n_rows;
		coord = [this, row_offset](unsigned j) { return m_mesh.yj(j + row_offset); };
		if (fixed != 0)
			fixed = n_cols - 1;
	}

	/* Case Neumann Boundary */
	if (m_bc_types.at(boundary) == 'n') {
		/* Initialize the vector that stores the boundary term values */
		m_nbc_matrix[boundary] = std::vector<double>(n, 0.);
		/* Assign values */
		for (; varying < n; ++varying)
			m_nbc_matrix[boundary][varying] = m_bc[boundary](coord(varying));
	}
	else {
		/* Make sure varying and fixed coordinates are in the right place when indexing m_solution */
		unsigned& x = (boundary % 2 == 0) ? varying : fixed;
		unsigned& y = (boundary % 2 == 0) ? fixed : varying;
		/* Assign boundary values */
		for (; varying < n; ++varying)
			m_solution[(y + not_rank0) * n_cols + x] = m_bc[boundary](coord(varying));
	}

}



/* Initialize m_solution in Parallel Cases */
void solver::parallel_init() const {

	/* Get mesh sizes */
	unsigned n_cols = m_mesh.nx(), n_rows = m_mesh.ny();

	/* Master Rank */
	if (rank == 0) {

		/* Size of solution = size of mesh so that it can then be collected more easily */
		m_solution = std::vector<double>(n_rows * n_cols, 0.);

		/* If the rows were perfectly evenly divided there's no need to find displacement and counts vectors */
		if (n_rows % size == 0)
			return;

		/* Initialize displacement and count vector if rows weren't evenly distributed */
		displacements.push_back(0);
		for (unsigned i = 0; i < size; ++i) {
			unsigned ith_rows = (i < (n_rows% size)) ? (n_rows / size) + 1 : (n_rows / size);
			counts_recv.push_back(ith_rows * n_cols);
			displacements.push_back(displacements.back() + counts_recv.back());
		}

	}
	else {
		/* Get amount of rows per rank + 1 for the extra bottom row */
		unsigned local_rows = (static_cast<unsigned>(rank) < (n_rows % size)) ? (n_rows / size) + 2 : (n_rows / size) + 1;
		/* Ranks other than size - 1 also have an extra row at the top */
		if (rank < size - 1)
			local_rows += 1;
		/* Initialize m_solution to be the correct size */
		m_solution = std::vector<double>(local_rows * n_cols, 0.);
		/* Free up the space that was occupied by default initialization */
		m_solution.shrink_to_fit();
	}

}



/* "Assigns" Neumann Boundary Conditions, i.e. gets the correct iteration update rule and the correct iteration bounds */
void solver::assign_nbc(const std::vector<unsigned>& boundaries, const std::vector<double>& prev_it, std::vector<unsigned>& limits, unsigned n_rows,
						const std::function<double(double, double, double, double, unsigned, unsigned)>& base, std::function<double(unsigned, unsigned)>& update_rule) const {

	/* Get number of columns since that is always the same for every rank */
	unsigned n_cols = m_mesh.nx();
	/* Ranks other than 0 have an extra row at the bottom so we need to take it into account when indexing */
	unsigned not_rank0 = (rank > 0 ? 1 : 0);
	/* Get step length in the x and y directions */
	double hx = m_mesh.delta_x(), hy = m_mesh.delta_y();

	/* Update iteration bounds on Neumann boundaries */
	for (unsigned b : boundaries) 
		if (m_bc_types.at(b) == 'n') {
			int offset = (b > 1) ? -1 : 1;
			limits[b] += offset;
		}

	/* Get correct update_rule that checks when we are on a boundary and uses the ghost node method to get a valid update */
	update_rule = [&prev_it, n_rows, n_cols, not_rank0, hx, hy, &base, this](unsigned i, unsigned j) {
		std::vector<double> v(4);

		if (j == n_rows - 1)				// Neumann BC on side 0
			v[0] = prev_it[(j - 1) * n_cols + i] + 2 * hy * m_nbc_matrix[0][i];
		else 
			v[0] = prev_it[(j + 1) * n_cols + i];

		if (i == n_cols - 1)				// Neumann BC on side 1
			v[1] = prev_it[j * n_cols + i - 1] + 2 * hx * m_nbc_matrix[1][j - not_rank0];
		else
			v[1] = prev_it[j * n_cols + i + 1];

		if (j == 0)							// Neumann BC on side 2
			v[2] = prev_it[(j + 1) * n_cols + i] - 2 * hy * m_nbc_matrix[2][i];
		else 
			v[2] = prev_it[(j - 1) * n_cols + i];
		
		if (i == 0)							// Neumann BC on side 3
			v[3] = prev_it[j * n_cols + i + 1] - 2 * hx * m_nbc_matrix[3][j - not_rank0];
		else
			v[3] = prev_it[j * n_cols + i - 1];

		return base(v[0], v[1], v[2], v[3], i, j);
	};

}



/* Sequential Solver */
void solver::sequential_solve(unsigned max_it, double tol) const {

	/* Ranks other than 0 won't have the correct m_solution size due to parallel_init() */
	if (rank > 0) {
		std::cerr << "Cannot call sequential_solve with a non-master rank" << std::endl;
		return;
	}

	/* Get mesh sizes */
	const unsigned n_cols = m_mesh.nx(), n_rows = m_mesh.ny();
	/* Variable to store the previous iteration */
	std::vector<double> prev_it = m_solution;
	/* Variable for the convergence check */
	double convergence_diff = 0.;
	/* Iteration Bounds for the base (Dirichlet) case */
	std::vector<unsigned> limits{n_rows - 2, n_cols - 2, 1, 1};
	/* Step length in the x and y directions */
	const double hx = m_mesh.delta_x(), hy = m_mesh.delta_y();
	/* Square of the tolerance to avoid calculating a square root each iteration */
	const double tol2 = tol * tol;

	/* Base of the update rule */
	std::function<double(double, double, double, double, unsigned, unsigned)> base;
	if (m_mesh.uniform())	// Uniform Mesh case
		base = [hx, n_cols, this](double v0, double v1, double v2, double v3, unsigned i, unsigned j) {return 0.25 * (v0 + v1 + v2 + v3 + hx * hx * m_force_matrix[j * n_cols + i]); };
	else					// Non Uniform Mesh case
		base = [hx, hy, n_cols, this](double v0, double v1, double v2, double v3, unsigned i, unsigned j) {return 0.5 * (hx * hx * (v0 + v2) + hy * hy * (v1 + v3) + hx * hx * hy * hy * m_force_matrix[j * n_cols + i]) / (hx * hx + hy * hy); };

	/* Actual update rule */
	std::function<double(unsigned, unsigned)> update_rule;

	if (m_bc_types.find('n') != std::string::npos)	// Case Neumann: Update Iteration bounds and get an update rule that checks when you are on a boundary
		assign_nbc(std::vector<unsigned>({ 0, 1, 2, 3 }), prev_it, limits, n_rows, base, update_rule);
	else											// Case Dirichlet: No need to check were we are since we can't be on a boundary
		update_rule = [&prev_it, n_cols, &base](unsigned i, unsigned j) {double v0 = prev_it[(j + 1) * n_cols + i], v1 = prev_it[j * n_cols + i + 1], v2 = prev_it[(j - 1) * n_cols + i], v3 = prev_it[j * n_cols + i - 1];
																		 return base(v0, v1, v2, v3, i, j); };

	/* Jacobi Iterations */
	for (unsigned k = 0; k < max_it; ++k) {
		for (unsigned j = limits[2]; j <= limits[0]; ++j) {
			for (unsigned i = limits[3]; i <= limits[1]; ++i) {
				m_solution[j * n_cols + i] = update_rule(i, j);							// Update each solution on each node
				double diff = m_solution[j * n_cols + i] - prev_it[j * n_cols + i];		// Get difference from previous iteration
				convergence_diff += diff * diff;										// Sum all these differences
			}
		}
		if (convergence_diff < tol2) {													// Convergence check
			converged = true;
			m_iterations = k + 1;
			break;
		}
		/* Reset for next iteration and put mark solution as previous */
		convergence_diff = 0.;
		std::swap(m_solution, prev_it);
	}

	/* If the solver hasn't converged it means we used up all iterations */
	m_iterations = (!converged) ? max_it : m_iterations;

}



/* Parallel Solver */
void solver::parallel_solve(unsigned max_it, double tol) const {

	/* Cannot call parallel solver if there are no other processes to communicate with */
	if (size < 2) {
		std::cerr << "Cannot call parallel_solve when size is less than 2" << std::endl;
		return;
	}

	/* Get mesh sizes */
	const unsigned n_cols = m_mesh.nx(), n_rows = m_mesh.ny();
	/* Get rows for each rank */
	unsigned local_rows = static_cast<unsigned>(rank) < (n_rows % size) ? (n_rows / size) + 1 : (n_rows / size);
	if (rank != 0 && rank != size - 1) 
		local_rows += 1; // +1 is needed because the row with index 0 isn't property of this rank AND the last row that's still property of this rank isn't actually a boundary one

	/* Variable to store the previous iteration */
	std::vector<double> prev_it = m_solution;
	/* Variable for the convergence check */
	double convergence_diff = 0.;
	/* Iteration Bounds for the base (Dirichlet) case */
	std::vector<unsigned> limits{ local_rows - 1, n_cols - 2, 1, 1 };
	/* Step length in the x and y directions */
	const double hx = m_mesh.delta_x(), hy = m_mesh.delta_y();
	/* Square of the tolerance to avoid calculating a square root each iteration */
	const double tol2 = tol * tol;

	/* Base of the update rule */
	std::function<double(double, double, double, double, unsigned, unsigned)> base;
	if (m_mesh.uniform())	// Uniform Mesh case
		base = [hx, n_cols, this](double v0, double v1, double v2, double v3, unsigned i, unsigned j) {return 0.25 * (v0 + v1 + v2 + v3 + hx * hx * m_force_matrix[j * n_cols + i]); };
	else 					// Non Uniform Mesh case
		base = [hx, hy, n_cols, this](double v0, double v1, double v2, double v3, unsigned i, unsigned j) {return 0.5 * (hx * hx * (v0 + v2) + hy * hy * (v1 + v3) + hx * hx * hy * hy * m_force_matrix[j * n_cols + i]) / (hx * hx + hy * hy); };

	/* Actual update rule */
	std::function<double(unsigned, unsigned)> update_rule;

	/* Neumann case, each rank should be treated separately */
	if (rank == 0 && (m_bc_types.at(1) == 'n' || m_bc_types.at(2) == 'n' || m_bc_types.at(3) == 'n'))				// Rank 0, only boundaries 1, 2 and 3 matter
			assign_nbc(std::vector<unsigned>({ 1, 2, 3 }), prev_it, limits, local_rows + 1, base, update_rule);
	else if (rank == size - 1 && (m_bc_types.at(0) == 'n' || m_bc_types.at(1) == 'n' || m_bc_types.at(3) == 'n'))	// Rank size - 1, only boundaries 0, 1 and 3 matter
			assign_nbc(std::vector<unsigned>({ 0, 1, 3 }), prev_it, limits, local_rows + 1, base, update_rule);
	else if (rank != 0 && rank != size - 1 && (m_bc_types.at(1) == 'n' || m_bc_types.at(3) == 'n'))					// Other Ranks, only boundaries 1 and 3 matter
			assign_nbc(std::vector<unsigned>({ 1, 3 }), prev_it, limits, local_rows + 1, base, update_rule);
	/* Dirichlet case */
	else												
		update_rule = [&prev_it, n_cols, &base](unsigned i, unsigned j) {double v0 = prev_it[(j + 1) * n_cols + i], v1 = prev_it[j * n_cols + i + 1], v2 = prev_it[(j - 1) * n_cols + i], v3 = prev_it[j * n_cols + i - 1];
																		 return base(v0, v1, v2, v3, i, j); };

	/* Jacobi Iterations */
	for (unsigned k = 0; k < max_it; ++k) {
		#pragma omp parallel for num_threads(4/size == 0 ? 1 : 4/size) reduction(+:convergence_diff)	// My machine gets worse performance from using threads differently from this
		for (unsigned j = limits[2]; j <= limits[0]; ++j) {
			for (unsigned i = limits[3]; i <= limits[1]; ++i) {
				m_solution[j * n_cols + i] = update_rule(i, j);											// Update each solution on each node
				double diff = m_solution[j * n_cols + i] - prev_it[j * n_cols + i];						// Get difference from previous iteration
				convergence_diff += diff * diff;														// Sum all these differences
			}
		}

		MPI_Allreduce(MPI_IN_PLACE, &convergence_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);			// Get global difference by having all ranks sum their local ones

		if (convergence_diff < tol2) {																	// Convergence check
			converged = true;
			m_iterations = k;
			break;
		}

		/* Communicate rows that are on the boundary between two ranks */
		if (rank > 0) {
			MPI_Recv(m_solution.data(), n_cols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(m_solution.data() + n_cols, n_cols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		}
		if (rank < size - 1) {
			MPI_Send(m_solution.data() + n_cols * (local_rows - 1), n_cols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(m_solution.data() + n_cols * local_rows, n_cols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		/* Reset for next iteration and put mark solution as previous */
		convergence_diff = 0.;
		std::swap(m_solution, prev_it);

	}


	/* Gather all local solutions into rank 0 */
	if (n_rows % size != 0) {	// Case with non-even distribution
		if (rank == 0) 
			MPI_Gatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, m_solution.data(), counts_recv.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		else if(rank == size - 1)
			MPI_Gatherv(m_solution.data() + n_cols, m_solution.size() - n_cols, MPI_DOUBLE, nullptr, counts_recv.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		else
			MPI_Gatherv(m_solution.data() + n_cols, m_solution.size() - 2*n_cols, MPI_DOUBLE, nullptr, counts_recv.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}	
	else {						// Case with even distribution
		if (rank == 0)
			MPI_Gather(MPI_IN_PLACE, local_rows * n_cols, MPI_DOUBLE, m_solution.data(), (n_rows / size)*n_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		else if (rank == size - 1)
			MPI_Gather(m_solution.data() + n_cols, m_solution.size() - n_cols, MPI_DOUBLE, nullptr, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		else
			MPI_Gather(m_solution.data() + n_cols, m_solution.size() - 2*n_cols, MPI_DOUBLE, nullptr, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	/* If the solver hasn't converged it means we used up all iterations */
	m_iterations = (!converged) ? max_it : m_iterations;

}



/* L2 Error */
double solver::L2_error(const Func& exact_solution) const {

	unsigned n_cols = m_mesh.nx(), n_rows = m_mesh.ny();

	double res = 0.;

	/* Calculate Error */
	for (unsigned i = 0; i < n_rows; ++i)
		for (unsigned j = 0; j < n_cols; ++j) {
			auto diff = (m_solution[i * n_cols + j] - exact_solution(m_mesh.xi(j), m_mesh.yj(i)));
			res += diff * diff;
		}

	return std::sqrt(m_mesh.delta_x() * m_mesh.delta_y() * res);

}




/* Write to VTK */
int solver::to_VTK(const std::string& filename) const {
	
	/* Write mesh to .vtk */
	int failure = m_mesh.to_VTK(filename);
	if (failure) {
		return 1;
	}

	/* Open the same file in append mode to add the solution */
	std::ofstream vtkFile(filename, std::ios_base::app);

	/* Check if the file opened */
	if (!vtkFile.is_open()) {
		std::cerr << "Error: could not open file " << filename << std::endl;
		return 1;
	}

	/* Write to file */
	vtkFile << "SCALARS scalars double\n";               // Description of the scalar field
	vtkFile << "LOOKUP_TABLE default\n";                 // Color table

	// Write vector field data
	for (const auto& u_ij : m_solution) {
		vtkFile << u_ij << "\n";
	}

	vtkFile.close();
	return 0;

}

}