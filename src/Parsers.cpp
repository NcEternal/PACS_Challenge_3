#include "Parsers.hpp"
#include "muParser.h"

namespace parsing {

/* ---------------------------------- */
/* ---------- Force Class ----------- */
/* ---------------------------------- */

/* Default Constructor */
Force::Force() {
	this->m_function = new mu::Parser;
	this->m_function->DefineVar("x", &m_x);
	this->m_function->DefineVar("y", &m_y);
}

/* Constructor given the function */
Force::Force(const std::string& expr) : Force() {
	m_expression = expr;
	this->m_function->SetExpr(expr);
}

/* Copy Constructor */
Force::Force(const Force& f) : m_function(new mu::Parser), m_x(f.m_x), m_y(f.m_y), m_expression(f.m_expression) {
	this->m_function->SetExpr(m_expression);
	this->m_function->DefineVar("x", &m_x);
	this->m_function->DefineVar("y", &m_y);
}

/* Assignment Operator */
Force& Force::operator=(const Force& f) {
	if (this != &f) {
		this->m_function->ClearVar();
		this->m_expression = f.m_expression;
		this->m_x = f.m_x;
		this->m_y = f.m_y;
		this->m_function->SetExpr(m_expression);
		this->m_function->DefineVar("x", &m_x);
		this->m_function->DefineVar("y", &m_y);
	}
	return *this;
}

/* Setter */
void Force::set_expression(const std::string& new_expr) {
	m_expression = new_expr;
	this->m_function->SetExpr(new_expr);
}

/* Like-a-function evaluation */
double Force::operator() (double x, double y) const {
	this->m_x = x;
	this->m_y = y;
	return this->m_function->Eval();
}

/* Destructor */
Force::~Force() {
	this->m_function->ClearVar();
	delete m_function;
}



/* ---------------------------------- */
/* ---- Boundary_Condition Class ---- */
/* ---------------------------------- */

/* Default Constructor */
Boundary_Condition::Boundary_Condition() {
	this->m_function = new mu::Parser;
	this->m_function->DefineVar("x", &m_v);
}

/* Constructor given the function */
Boundary_Condition::Boundary_Condition(const std::string& expr) : Boundary_Condition() {
	m_expression = expr;
	this->m_function->SetExpr(expr);
}

/* Constructor given the function and its variable */
Boundary_Condition::Boundary_Condition(const std::string& expr, char var_name) : m_expression(expr), m_var_name(var_name) {
	this->m_function = new mu::Parser;
	this->m_function->SetExpr(expr);
	this->m_function->DefineVar(std::string({ m_var_name }), &m_v);
}

/* Copy Constructor */
Boundary_Condition::Boundary_Condition(const Boundary_Condition& bc) : m_function(new mu::Parser), m_v(bc.m_v), m_expression(bc.m_expression), m_var_name(bc.m_var_name) {
	this->m_function->SetExpr(m_expression);
	this->m_function->DefineVar(std::string({ m_var_name }), &m_v);
}

/* Assignment Operator */
Boundary_Condition& Boundary_Condition::operator=(const Boundary_Condition& bc) {
	if (this != &bc) {
		this->m_function->ClearVar();
		this->m_expression = bc.m_expression;
		this->m_v = bc.m_v;
		this->m_var_name = bc.m_var_name;
		this->m_function->SetExpr(m_expression);
		this->m_function->DefineVar(std::string({ m_var_name }), &m_v);
	}
	return *this;
}

/* Setter given new function with the same variable as before */
void Boundary_Condition::set_expression(const std::string& new_expr) {
	m_expression = new_expr;
	this->m_function->SetExpr(new_expr);
}

/* Setter given new function with a different variable from before */
void Boundary_Condition::set_expression(const std::string& new_expr, char new_var_name) {
	m_var_name = new_var_name;
	m_expression = new_expr;
	this->m_function->SetExpr(new_expr);
	this->m_function->ClearVar();
	this->m_function->DefineVar(std::string({ m_var_name }), &m_v);
}

/* Like-a-function evaluation */
double Boundary_Condition::operator() (double var) const {
	this->m_v = var;
	return this->m_function->Eval();
}

/* Destructor */
Boundary_Condition::~Boundary_Condition() {
	this->m_function->ClearVar();
	delete m_function;
}

}