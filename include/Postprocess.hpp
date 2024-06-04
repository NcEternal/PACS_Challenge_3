#pragma once
#ifndef HH_Postprocess_HH
#define HH_Postprocess_HH

#include "Mesh.hpp"
#include <functional>
#include <vector>

namespace Errors {

/* L2 Error between numerical_solution and exact_solution over mesh */
double L2_error(const std::vector<double> numerical_solution, const std::function<double(double, double)>& exact_solution, const Mesh& mesh);

}



/* Function to write solution to .vtk file */
void to_VTK(const std::vector<double>& solution, const Mesh& mesh, const std::string& filename);

#endif