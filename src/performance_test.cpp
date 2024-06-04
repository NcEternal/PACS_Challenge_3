#include "Solver.hpp"

#include <mpi.h>	// parallelism

#include <sstream>	// create file name
#include <fstream>	// write to file
#include <iomanip>	// setfill, setw

#include <numbers>	// pi
#include <cmath>	// sin

#include <chrono>	// timing

#ifndef REPS
#define REPS 1
#endif

/* Performance test over a default problem */

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	auto pi = std::numbers::pi;
	auto f = [pi](double x, double y) {return 8 * pi * pi * std::sin(2 * pi * x) * std::sin(2 * pi * y); };
	auto u_ex = [pi](double x, double y) {return std::sin(2 * pi * x) * std::sin(2 * pi * y); };

	std::ostringstream oss;
	oss << "../test/data/stats_" << size << "procs.csv";

	std::ofstream output_file(oss.str(), std::ofstream::out);
	if (rank == 0) {
		output_file << "n,calculation_time_(ms),iterations,L2_error,has_converged\n";
	}

	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();
	auto time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);

	for (unsigned k = 4; k < 9; ++k) {

		unsigned n = 1 << k;
		Mesh mesh(0, 1, 0, 1, n);
		double avg = 0.;

		for (unsigned r = 0; r < REPS - 1; ++r) {

			laplace::solver ls(mesh, f, rank, size);

			if (size < 2) {
				t1 = std::chrono::high_resolution_clock::now();
				ls.sequential_solve(1e7, 1e-6);
				t2 = std::chrono::high_resolution_clock::now();
				time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
			}
			else {
				t1 = std::chrono::high_resolution_clock::now();
				ls.parallel_solve(1e7, 1e-6);
				t2 = std::chrono::high_resolution_clock::now();
				time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
			}

			double diff = time_taken.count() / 1000000. - avg;
			avg += diff / (r + 1);

		}

		laplace::solver ls(mesh, f, rank, size);

		if (size < 2) {
			t1 = std::chrono::high_resolution_clock::now();
			ls.sequential_solve(1e7, 1e-6);
			t2 = std::chrono::high_resolution_clock::now();
			time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
		}
		else {
			t1 = std::chrono::high_resolution_clock::now();
			ls.parallel_solve(1e7, 1e-6);
			t2 = std::chrono::high_resolution_clock::now();
			time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
		}

		double diff = time_taken.count() / 1000000. - avg;
		avg += diff / REPS;

		if (rank == 0) {
			output_file << n << "," << avg << "," << ls.iterations() << "," << ls.L2_error(u_ex) << "," << ls.has_converged() << "\n";

			std::ostringstream ossvtk;
			ossvtk << "../test/data/solution_" << size << "procs_" << std::setfill('0') << std::setw(5) << n << ".vtk";
			ls.to_VTK(ossvtk.str());
		}
		if (size > 1)
			MPI_Barrier(MPI_COMM_WORLD);

	}

	output_file.close();
	MPI_Finalize();

	return 0;

}