# Challenge 3 A Matrix-free Parallel Solver for the Laplace Equation #

Solution to the $3^{rd}$ PACS challenge. Implementation of a Laplace Equation solver.

The make file produces two executables: `main` and `performance_test`. <br>
`main` is the main way to use the solver. It can be run with from the command line using `mpiexec -np n ./main`, where n is the number of
processes you want to use. Upon launching, it will ask for the name of a file (extension included) in the `/src/parameters/`
folder. Pressing enter without providing such name will use the default problem data. The solution to the equation will be
written in *.vtk* format in the `/src/output/` folder in a file named `out_{name of the input file without extension}.vtk` <br>
`performance_test` is used to test the speed and precision of the class on the default problem
for a grid with $n = 2^k, k = 4, ..., 8$ points per side. `performance_test` is best use with the script
found in the `./test/` directory.



## The Default Problem ##

The default problem in both `main` and `performance_test` is the following:

$$
	\begin{cases}
		- \Delta u = f(x), & \quad \text{in $\Omega = (0,1) \times (0,1)$,} \\
		u=0, & \quad \text{on $\{x=0\}$,}\\
		u=0, & \quad \text{on $\{x=1\}$,}\\
		u=0, & \quad \text{on $\{y=0\}$,}\\
		u=0, & \quad \text{on $\{y=1\}$}\\
	\end{cases}
$$



## The `parameters` Struct ##

Inside the namespace `laplace` is defined the struct `parameters` which collects all parameters
relative to the problem, initialized by default to those of the above mentioned problem. The struct
can be found in the `Parameters.hpp` file together with an explanation of what each parameter does



## The `Mesh` Class ##

Inside the file `Mesh.hpp` is defined the `Mesh` calss, representing a Cartesian decomposition of a 
rectangle. It comes with two constructors:

* `Mesh(double xlb, double xub, double ylb, double yub, unsigned nx_, unsigned ny_)`;

* `Mesh(double xlb, double xub, double ylb, double yub, unsigned n_)`;

where `xlb` and `ylb` are the lower bounds of the intervals defining the rectangle, `xub` and `yub` the
upper bounds, and `nx_` and `ny_` represent how many points each interval should be divided in. The
constructor with only `n_` imposes `nx_` = `ny_` = `n_`. <br>

The class provides the following methods:

* `double xi(unsigned i) const`: returns the point $x_i$ of the decomposition (where $x_0$ = `xlb` and $x_{nx-1}$ = `xub`);

* `double yj(unsigned j) const`: returns the point $y_j$ of the decomposition (where $y_0$ = `ylb` and $y_{ny-1}$ = `yub`);

* `unsigned nx() const`: returns `nx_`;

* `unsigned ny() const`: returns `ny_`;

* `double delta_x() const`: returns the distance between consecutive points along the x axis;

* `double delta_y() const`: returns the distance between consecutive points along the y axis;

* `bool uniform() const`: returns true if `delta_x()` = `delta_y()`;

* `int to_VTK(const std::string& filename) const`: prints the information about the mesh in *.vtk* format on the
file `filename`. Returns 1 if it couldn't open the file, 0 otherwise. 



## The `solver` Class ##

Inside the `Solver.hpp` file, in the `laplace` namespace is defined the `solver` class. This is the class used to 
numerically solve the Laplace equation using the Jacobi iteration method, either sequentially or parallely. 
Going forward we'll use the following aliases: `Func = std::function<double(double, double)>` and `BC_Func = std::function<double(double)>`;

The class comes with the following constructors:

1. `solver(const Func& force)`: initializes the Laplace equation with forcing term `force` over the default domain
with the default boundary conditions and `nx` = `ny` = 32;

2. `solver(unsigned nx, unsigned ny, const Func& force)`: same as (1), but where `nx` and `ny` can be chosen individually;

3. `solver(unsigned n, const Func& force)`: same as (1), but where `nx` = `ny` = `n`;

4. `solver(const Mesh& mesh, const Func& force)`: initializes the Laplace equation with forcing term `force` over the domain described
by `mesh` with default boundary conditions;

5. `solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned nx, unsigned ny, const Func& force)`: same as (4), but the mesh is
constructed automatically with the given parameters;

6. `solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned n, const Func& force)`: same as (5);

7. `solver(const Mesh& mesh, const Func& force, const BC_Func& bc0, const BC_Func& bc1, const BC_Func& bc2, const BC_Func& bc3, const std::string& bc_types)`:
initialized the Laplace equation with forcing term `force` over the domain described by `mesh` and boundary conditions described by
`bc0`, `bc1`, `bc2`, `bc3` and `bc_types`. More in depth information on how the boundary conditions are defined can be found in both the `Parameters.hpp` file
and the `sample.txt` file in the `/src/parameters/` folder;

8. `solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned nx, unsigned ny, const Func& force, const BC_Func& bc0, const BC_Func& bc1, const BC_Func& bc2, const BC_Func& bc3, const std::string& bc_types)`:
same as (7), but the mesh is constructed automatically with the given parameters;

9. `solver(double x_lb, double x_ub, double y_lb, double y_ub, unsigned n,  const Func& force, const BC_Func& bc0, const BC_Func& bc1, const BC_Func& bc2, const BC_Func& bc3, const std::string& bc_types)`:
same as (8)

All these constructors come with a parallel version that takes two extra arguments at the end: `int rank_` and `int size_` which represent the rank
of the process calling the constructor and the number of processes respectively.<br>

The class provides the following methods:

* `void sequential_solve(unsigned max_it, double tol) const`: sequentially solves the Laplace equation. `max_it` is the maximum number of iterations
allowed while `tol` is such that the algorithm is considered to have converged if the euclidean norm of the difference between two iterations is less than
`tol`. `sequential_solve` can only be called by objects that were initialized with no given `rank` or with `rank` = 0;

* `void parallel_solve(unsigned max_it, double tol) const`: parallely solves the Laplace equation. `max_it` and `tol` have the same roles as before.
`parallel_solve` can only be called if the object was initialized with `size` > 1;

* `std::vector<double> solution() const`: returns the solution after a call to one of the solving methods. The solution is such that the value stored at
the index `j*nx + i` corresponds to the value of the solution at the point $(x_i, y_j)$;

* `unsigned iterations() const`: returns the amount of iterations taken by the last call to one of the solving methods;

* `bool has_converged() const`: returns `true` if the last call to a solving method converged, `false` otherwise;

* `Mesh mesh() const`: returns the mesh used by the object;

* `Func force() const`: returns the force used by the object;

* `std::vector<BC_Func> boundary_conditions() const`: returns the 4 boundary conditions used by the object;

* `BC_Func boundary_condition(unsigned i) const`: returns the i-th boundary condition used by the object;

* `std::string boundary_types() const`: returns the string describing the boundary types used by the object;

* `void set_force(const Func& f)`: changes the forcing term to be equal to `f`. Resets convergence status, number of iterations and solution;

* `void set_boundary_conditions(const std::vector<BC_Func>& bcs)`: changes the boundary conditions to be equal to `bcs`. Resets
convergence status, number of iterations and solution;

* `void set_boundary_condition(unsigned i, const BC_Func& bc)`: changes the i-th boundary condition to be equal to `bc`. Resets
convergence status, number of iterations and solution;

* `void set_boundary_types(const std::string& bc_types)`: changes the string describing the boundary types to be equal to `bc_types`. 
Resets convergence status, number of iterations and solution;

* `double L2_error(const Func& exact_solution) const`: returns the $L^2$ error of the current solution given the exact one;

* `int to_VTK(const std::string& filename) const`: prints everything about the solution (mesh data included) in `filename` in *.vtk* format.
Returns 1 if it couldn't open the file, 0 otherwise. 



## Parsers ##

The `solver` class cannot accept `std::string` in place of one of the forces, thus parsers are provided in the `Parsers.hpp` file
under the namespace `parsing` to make using the `solver` class easier. There are two types of parsers:

* `Force`: represents a 2 variable function and can thus be used for both the forcing term and the exact solution. It comes with
the following:
	- `Force()`: the default constructor that initializes the object without storing any function in it;
	
	- `Force(const std::string& expr)`: constructor where `expr` is the function to be parsed;

	- `Force(const Force& f)`: copy constructor;

	- `Force& operator=(const Force& f)`: assignment operator;
		
	- `void set_expression(const std::string& new_expr)`: changes the function to be parsed to match `new_expr`;
		
	- `double operator() (double x, double y) const`: evaluates the parsed function at $(x, y)$;

* `Boundary_Condition`: represents a 1 variable function and can thus be used for the boundary terms. It comes with ù
the following:
	- `Boundary_Condition()`: the default constructor that initializes the object without storing any function in it;

	- `Boundary_Condition(const std::string& expr)`: constructor where `expr` is a function in **x** that needs to be parsed;

	- `Boundary_Condition(const std::string& expr, char var_name)`: constructor where `expr` is a function in
		`var_name` that needs to be parsed;

	- `Boundary_Condition(const Boundary_Condition& bc)`: the copy constructor;

	- `Boundary_Condition& operator=(const Boundary_Condition& bc)`: the assignment operator;

	- `void set_expression(const std::string& new_expr)`: changes the function to be parsed to match `new_expr` assuming
	it is in the same variable as the one before. If default constructed, the variable is assumed to be **x**;

	- `void set_expression(const std::string& new_expr, char new_var_name)`: changes the function to be parsed to match `new_expr` 
	assuming it is in the variable `new_var_name`;

	- `double operator() (double var) const`: evaluates the parsed function when its variable equals `var`.

Both parser make use of the **muParser** library. As such the strings passed to them should be **muParser** readable.



## Input and Output ##

In the file `Input_Output`, under the namespace `laplace`, are defined the following functions:

* `parameters read(const std::string& filename)`: reads file `filename` and returns a `parameter` struct matching
the information in the file. For more in depth information on the file format look at the `sample.txt` file in the
`/src/parameters/` folder; 

* `parameters read(const std::string& filename, int rank)`: an overload of the previous function that lets the caller
read the file only if `rank == 0`;

* `void print_data(const parameters& params)`: outputs the problem data to `std::cout`;

* `void print_data(const parameters& params, int rank)`: an overload of the previous function that lets the caller
output the data only if `rank == 0`;



## Communication Between Processes ##

In the file `Parameter_Communication.hpp`, under the namespace `laplace`, are defined the following functions:

* `void communicate(parameters& p, int rank, int size)`: if `rank == 0` it communicates the data contained in `p`
to the other processes (if there are any). if `rank != 0` it receives parameter data from 0 and stores it in `p`;

* `void communicate(parameters& p)`: an overload of the previous function that gets the values of `rank` and `size`
on its own.



## Utilities ##

In the file `Postprocess.hpp` are defined the following functions:

* `double L2_error(const std::vector<double> numerical_solution, const std::function<double(double, double)>& exact_solution, const Mesh& mesh)`:
inside of the namespace `Errors`, it provides an out of class alternative to `solver::L2_error`;

* `void to_VTK(const std::vector<double>& solution, const Mesh& mesh, const std::string& filename)`: it provides an out of class alternative to
`solver::to_VTK`.



## Adding New Equations to be Solved ##

To add new problems, read the `sample.txt` file in the `/src/parameters/` folder, create a file following its format
and add it to the folder. In the same folder there are multiple other problems already
implemented. `messy.txt` provides an exampe of how bad a readable file can get.