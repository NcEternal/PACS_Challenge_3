#include "Postprocess.hpp"
#include <cmath>
#include <fstream>
#include <iostream>

namespace Errors {

double L2_error(const std::vector<double> numerical_solution, const std::function<double(double, double)>& exact_solution, const Mesh& mesh) {

	/* Get Mesh sizes */
	unsigned n_cols = mesh.nx(), n_rows = mesh.ny();

	double res = 0.;

	/* Calculate Error */
	for (unsigned i = 0; i < n_rows; ++i)
		for (unsigned j = 0; j < n_cols; ++j) {
			auto diff = (numerical_solution[i * n_cols + j] - exact_solution(mesh.xi(j), mesh.yj(i)));
			res += diff * diff;
		}

	return std::sqrt(mesh.delta_x() * mesh.delta_y() * res);

}

}



/* Function to write solution to .vtk file */
void to_VTK(const std::vector<double>& solution, const Mesh& mesh, const std::string& filename) {

	/* Write Mesh to VTK */
	int failure = mesh.to_VTK(filename);
	if (failure) {
		return;
	}

	std::ofstream vtkFile(filename, std::ios_base::app);

	if (!vtkFile.is_open()) {
		std::cerr << "Error: could not open file " << filename << std::endl;
		return;
	}

	vtkFile << "SCALARS scalars double\n";               // description of the scalar field
	vtkFile << "LOOKUP_TABLE default\n";                 // color table

	/* Write vector field data */
	for (const auto& u_ij : solution) {
		vtkFile << u_ij << "\n";
	}

	vtkFile.close();

}