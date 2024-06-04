#pragma once
#ifndef HH_Parsers_HH
#define HH_Parsers_HH

#include <string>

/* Forward Declaration to avoid including muParser.hpp in the header file */
namespace mu {
class Parser;
}

namespace parsing {

/* Class to represent a 2 variable function (usually the force of a Laplace Problem, but also its solution) */

class Force {

private:
	mu::Parser* m_function;
	mutable double m_x;
	mutable double m_y;
	std::string m_expression;

public:
	/* Constructors */
	Force();
	Force(const std::string& expr);
	
	/* Copy Constructor */
	Force(const Force& f);
	
	/* Assignment Operator */
	Force& operator=(const Force& f);

	/* Setter */
	void set_expression(const std::string& new_expr);

	/* Like-a-function evaluation */
	double operator() (double x, double y) const;

	/* Destructor */
	~Force();

};



/* Class to represent a 1 variable function (usually a boundary condition of a Laplace Problem) */

class Boundary_Condition {

private:
	mu::Parser* m_function;
	mutable double m_v;
	std::string m_expression;
	char m_var_name = 'x';

public:
	/* Constructors */
	Boundary_Condition();
	Boundary_Condition(const std::string& expr);
	Boundary_Condition(const std::string& expr, char var_name);

	/* Copy Constructor */
	Boundary_Condition(const Boundary_Condition& bc);

	/* Assignment Operator */
	Boundary_Condition& operator=(const Boundary_Condition& bc);

	/* Setters */
	void set_expression(const std::string& new_expr);
	void set_expression(const std::string& new_expr, char new_var_name);

	/* Like-a-function evaluation */
	double operator() (double var) const;

	/* Destructor */
	~Boundary_Condition();

};

}

#endif