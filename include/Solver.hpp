#pragma once
#ifndef HH_Solver_HH
#define HH_Solver_HH

#include "Mesh.hpp"
#include <functional>	
#include <vector>		

namespace laplace{

/* Class to numerically solve the Laplace Problem with any force/boundary condition of types Dirichlet/Neumann over any rectangle */

class solver {

	using Func = std::function<double(double, double)>;
	using BC_Func = std::function<double(double)>;

private:
	/* Problem Data */
	Mesh m_mesh{ Mesh(0, 1, 0, 1, 32) };
	Func m_force;
	std::vector<BC_Func> m_bc{ ([](double x) {return 0.; }), ([](double y) {return 0.; }), ([](double x) {return 0.; }), ([](double y) {return 0.; }) };
	std::string m_bc_types{ "dddd" };
	mutable std::vector<double> m_force_matrix;
	mutable std::vector<std::vector<double>> m_nbc_matrix{ std::vector<double>(), std::vector<double>(), std::vector<double>(), std::vector<double>() };
	
	/* Solution */
	mutable std::vector<double> m_solution{ std::vector<double>(m_mesh.ny() * m_mesh.nx(), 0.) };
	
	/* Convergence */
	mutable unsigned m_iterations = 0;
	mutable bool converged = false;
	
	/* Parallel Variables */
	const unsigned rank = 0;
	const unsigned size = 1;
	mutable std::vector<int> counts_recv, displacements;	

	/* Note: counts_recv, displcements, m_force_matrix and m_nbc_matrix are mutable because they cannot be default initialized but we still want to use const solvers
	   On the same note, if m_solution, m_iterations or m_converged weren't mutable, const solvers would be useless from the start */



	/* Pre-calculates Force over the Mesh */
	void pre_calc_force() const;

	/* Assigns Boundary Conditions */
	void assign_bc() const;
	void assign_bc_helper(const unsigned boundary, const unsigned n_rows, const unsigned n_cols, const unsigned row_offset, const unsigned not_rank0) const;

	/* Initialize m_solution in Parallel Cases */
	void parallel_init() const;

	/* "Assigns" Neumann Boundary Conditions */
	void assign_nbc(const std::vector<unsigned>& boundaries, const std::vector<double>& prev_it, std::vector<unsigned>& limits, unsigned n_rows,
					const std::function<double(double, double, double, double, unsigned, unsigned)>& base, std::function<double(unsigned, unsigned)>& update_rule) const;



public:
	/* Sequential Constructors */
	solver(const Func& force) : m_force(force) { pre_calc_force(); }

	solver(const Mesh& mesh, const Func& force) : m_mesh(mesh), m_force(force) { pre_calc_force(); }
	
	solver(unsigned nx, unsigned ny, const Func& force) : solver(Mesh(0, 1, 0, 1, nx, ny), force) {}
	
	solver(unsigned n, const Func& force) : solver(Mesh(0, 1, 0, 1, n), force) {}
	
	solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned nx, unsigned ny, const Func& force) : solver(Mesh(x_lb, x_ub, y_lb, y_ub, nx, ny), force) {}
	
	solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned n, const Func& force) : solver(Mesh(x_lb, x_ub, y_lb, y_ub, n), force) {}
	
	solver(const Mesh& mesh, const Func& force, const BC_Func& bc0, const BC_Func& bc1, const BC_Func& bc2, const BC_Func& bc3, const std::string& bc_types) : 
		m_mesh(mesh), m_force(force), m_bc{ bc0, bc1, bc2, bc3 }, m_bc_types(bc_types) { assign_bc(); pre_calc_force(); }

	solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned nx, unsigned ny, const Func& force, const BC_Func& bc0, const BC_Func& bc1, const BC_Func& bc2, const BC_Func& bc3, const std::string& bc_types) :
		solver(Mesh(x_lb, x_ub, y_lb, y_ub, nx, ny), force, bc0, bc1, bc2, bc3, bc_types) {}

	solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned n, const Func& force, const BC_Func& bc0, const BC_Func& bc1, const BC_Func& bc2, const BC_Func& bc3, const std::string& bc_types) :
		solver(Mesh(x_lb, x_ub, y_lb, y_ub, n), force, bc0, bc1, bc2, bc3, bc_types) {}

	/* Parallel Constructors */
	solver(const Func& force, int rank_, int size_) : m_force(force), rank(rank_), size(size_) { parallel_init(); pre_calc_force(); }

	solver(const Mesh& mesh, const Func& force, int rank_, int size_) : m_mesh(mesh), m_force(force), rank(rank_), size(size_) { parallel_init(); pre_calc_force(); }

	solver(unsigned nx, unsigned ny, const Func& force, int rank_, int size_) : solver(Mesh(0, 1, 0, 1, nx, ny), force, rank_, size_) {}

	solver(unsigned n, const Func& force, int rank_, int size_) : solver(Mesh(0, 1, 0, 1, n), force, rank_, size_) {}

	solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned nx, unsigned ny, const Func& force, int rank_, int size_) : solver(Mesh(x_lb, x_ub, y_lb, y_ub, nx, ny), force, rank_, size_) {}

	solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned n, const Func& force, int rank_, int size_) : solver(Mesh(x_lb, x_ub, y_lb, y_ub, n), force, rank_, size_) {}

	solver(const Mesh& mesh, const Func& force, const BC_Func& bc0, const BC_Func& bc1, const BC_Func& bc2, const BC_Func& bc3, const std::string& bc_types, int rank_, int size_) :
		m_mesh(mesh), m_force(force), m_bc{ bc0, bc1, bc2, bc3 }, m_bc_types(bc_types), rank(rank_), size(size_) { parallel_init(); assign_bc(); pre_calc_force(); }

	solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned nx, unsigned ny, const Func& force, const BC_Func& bc0, const BC_Func& bc1, const BC_Func& bc2, const BC_Func& bc3,
				   const std::string& bc_types, int rank_, int size_) :
		solver(Mesh(x_lb, x_ub, y_lb, y_ub, nx, ny), force, bc0, bc1, bc2, bc3, bc_types, rank_, size_) {}

	solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned n, const Func& force, const BC_Func& bc0, const BC_Func& bc1, const BC_Func& bc2, const BC_Func& bc3,
				   const std::string& bc_types, int rank_, int size_) :
		solver(Mesh(x_lb, x_ub, y_lb, y_ub, n), force, bc0, bc1, bc2, bc3, bc_types, rank_, size_) {}



	/* Getters */
	std::vector<double> solution() const { return m_solution; }
	unsigned iterations() const { return m_iterations; }
	bool has_converged() const { return converged; }
	Mesh mesh() const { return m_mesh; }
	Func force() const { return m_force; }
	std::vector<BC_Func> boundary_conditions() const { return m_bc; }
	BC_Func boundary_condition(unsigned i) const { return m_bc[i % 4]; }
	std::string boundary_types() const { return m_bc_types; }

	/* Setters */	/* All of these reset the solution */
	void set_force(const Func& f) { m_force = f; converged = false; pre_calc_force(); m_solution = std::vector<double>(m_solution.size(), 0.); m_iterations = 0; }
	void set_boundary_conditions(const std::vector<BC_Func>& bcs) { m_bc = bcs; assign_bc(); converged = false; m_iterations = 0; }
	void set_boundary_condition(unsigned i, const BC_Func& bc) { m_bc[i % 4] = bc; assign_bc(); converged = false; m_iterations = 0; }
	void set_boundary_types(const std::string& bc_types) { m_bc_types = bc_types; assign_bc(); converged = false; m_iterations = 0; }

	/* Solvers */
	void sequential_solve(unsigned max_it, double tol) const;
	void parallel_solve(unsigned max_it, double tol) const;

	/* Errors */
	double L2_error(const Func& exact_solution) const;

	/* Write to VTK */
	int to_VTK(const std::string& filename) const;

};

}

#endif