#include "Input_Output.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace laplace {

	namespace helper {
	/* Helper function to turn the file into a parseable string */
	std::string file_to_string(const std::string& filename) {
		std::ifstream input_file(filename);

		if (!input_file.is_open()) {
			std::cerr << "Error: could not open file " << filename << std::endl;
			return std::string();
		}

		std::stringstream str_file;
		std::string line;
		std::string standardized_line;

		/* Read line by line */
		while (std::getline(input_file, line)) {
			/* Remove unwanted characters */
			std::erase(line, ' ');
			std::erase(line, '\t');
			std::erase(line, '\n');
			std::erase(line, '\r');
			/* Remove comments */
			standardized_line = line.substr(0, line.find('#'));
			/* If there was a field in the line, append '#' to delimit it */
			if (!standardized_line.empty())
				standardized_line.append("#");
			/* Add the line to the string representing the whole file */
			str_file << standardized_line;
		}

		return str_file.str();
	}
	}



/* Funtion to read the parameters of the problem from a file */
parameters read(const std::string& filename) {

	parameters default_params;

	/* No file to read ==> default data */
	if (filename.empty())
		return default_params;

	/* Turn the file into a string that can be more easily parsed */
	std::string p_str = helper::file_to_string(filename);
	/* File contained no data ==> default data*/
	if (p_str.empty())
		return default_params;

	/* Fields to search */
	const std::vector<std::string> search = { "x_lb=", "x_ub=", "y_lb=", "y_ub=", "nx=", "ny=", "force=", "bc_0=", "bc_1=", "bc_2=", "bc_3=", "bc_types=", "max_it=", "tol=", "u_ex="};
	/* Values on each field */
	std::vector<std::string> v_p_str(search.size());

	/* Get values of each field */
	for (unsigned i = 0; i < search.size(); ++i) {
		std::string target = search[i];
		auto loc = p_str.find(target);		// Find field
		if (loc != std::string::npos) {		// If found get its data 
			auto start = loc + target.size();
			v_p_str[i] = p_str.substr(start, p_str.find('#', loc) - start);		// file_to_string terminates each field with '#'
		}
	}

	/* Get actual values of the fields (unless they were left empty, in which case default values are used) */
	if (!v_p_str[0].empty())
		default_params.x_lb = std::stod(v_p_str[0]);
	if (!v_p_str[1].empty())
		default_params.x_ub = std::stod(v_p_str[1]);
	if (!v_p_str[2].empty())
		default_params.y_lb = std::stod(v_p_str[2]);
	if (!v_p_str[3].empty())
		default_params.y_ub = std::stod(v_p_str[3]);
	if (!v_p_str[4].empty())
		default_params.nx = static_cast<unsigned>(std::stoi(v_p_str[4]));
	if (!v_p_str[5].empty())
		default_params.ny = static_cast<unsigned>(std::stoi(v_p_str[5]));
	if (!v_p_str[6].empty())
		default_params.force = v_p_str[6];
	for (unsigned i = 0; i < 4; ++i) {
		if (!v_p_str[7 + i].empty())
			default_params.bcs[i] = v_p_str[7 + i];
	}
	if (!v_p_str[11].empty())
		default_params.bc_types = v_p_str[11];
	if (!v_p_str[12].empty())
		default_params.max_it = static_cast<unsigned>(std::stoi(v_p_str[12]));
	if (!v_p_str[13].empty())
		default_params.tol = std::stod(v_p_str[13]);
	if (!v_p_str[14].empty())
		default_params.exact_sol = v_p_str[14];

	return default_params;

}

/* Overload of the previous function that lets only rank 0 read the file */
parameters read(const std::string& filename, int rank) {
	parameters default_params;
	if (rank == 0)
		return read(filename);
	return default_params;
}



/* Function to print the parameters of the problem on the std::cout */
void print_data(const parameters& params) {
	std::cout << "Problem Data:\n";
	std::cout << "Domain: [" << params.x_lb << ", " << params.x_ub << "] x [" << params.y_lb << ", " << params.y_ub << "]\n";
	std::cout << "Refinements:\n\t- " << params.nx << " points on the x-axis;\n\t- " << params.ny << " points on the y-axis;\n";
	std::cout << "Forcing Term: " << params.force << "\n";
	std::cout << "Boundary Terms:\n\t- " << params.bcs[0] << " on y = " << params.y_ub << " (" << ((params.bc_types.at(0) == 'd') ? "Dirichlet" : "Neumann") << ");\n\t- " 
			  << params.bcs[1] << " on x = " << params.x_ub << " (" << ((params.bc_types.at(1) == 'd') ? "Dirichlet" : "Neumann") << ");\n\t- "
			  << params.bcs[2] << " on y = " << params.y_lb << " (" << ((params.bc_types.at(2) == 'd') ? "Dirichlet" : "Neumann") << ");\n\t- "
			  << params.bcs[3] << " on x = " << params.x_lb << " (" << ((params.bc_types.at(3) == 'd') ? "Dirichlet" : "Neumann") << ");\n";
	std::cout << "Maximum Iterations: " << params.max_it << "\n";
	std::cout << "Tolerance: " << params.tol << "\n" << std::endl;
}



/* Overload of the previous function that lets only rank 0 write to std::cout */
void print_data(const parameters& params, int rank) {
	if (rank == 0)
		print_data(params);
}

}