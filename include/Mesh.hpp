#pragma once
#ifndef HH_Mesh_HH
#define HH_Mesh_HH

#include <string>

/* Class to represent a Cartesian Decomposition of a rectangle */

class Mesh {
private:
	double x_lb, x_ub, y_lb, y_ub;
	unsigned m_nx, m_ny;
	double hx, hy;

public:
	/* Constructors */
	Mesh(double xlb, double xub, double ylb, double yub, unsigned nx_, unsigned ny_) : x_lb(xlb), x_ub(xub), y_lb(ylb), y_ub(yub), m_nx(nx_), m_ny(ny_) 
																					{ hx = (x_ub - x_lb) / (m_nx - 1), hy = (y_ub - y_lb) / (m_ny - 1); }
	Mesh(double xlb, double xub, double ylb, double yub, unsigned n_) : Mesh(xlb, xub, ylb, yub, n_, n_) {}

	/* Getters */
	double xi(unsigned i) const { return x_lb + i * hx; }
	double yj(unsigned j) const { return y_lb + j * hy; }
	unsigned nx() const { return m_nx; }
	unsigned ny() const { return m_ny; }
	double delta_x() const { return hx; }
	double delta_y() const { return hy; }
	
	/* Uniformity of the mesh */
	bool uniform() const { return hx == hy; }

	/* Write to VTK */
	int to_VTK(const std::string& filename) const;

};

#endif