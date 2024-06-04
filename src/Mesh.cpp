#include "Mesh.hpp"
#include <fstream>
#include <iostream>



/* Write Mesh to VTK */
int Mesh::to_VTK(const std::string& filename) const {

	std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return 1;
    }

    /* Write VTK header */
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Scalar Field Data\n";
    vtkFile << "ASCII\n";

    /* Write Mesh data */
    vtkFile << "DATASET STRUCTURED_POINTS\n";                               // format of the dataset
    vtkFile << "DIMENSIONS " << m_nx << " " << m_ny << " " << 1 << "\n";    // number of points in each direction
    vtkFile << "ORIGIN " << x_lb << " " << y_lb << " 0\n";                  // lower-left corner of the structured grid
    vtkFile << "SPACING" << " " << hx << " " << hy << " " << 1 << "\n";     // spacing between points in each direction
    vtkFile << "POINT_DATA " << m_nx * m_ny << "\n";                        // number of points

    vtkFile.close();

    return 0;

}
