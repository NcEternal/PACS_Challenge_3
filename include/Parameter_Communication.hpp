#pragma once
#ifndef HH_Parameter_Communication_HH
#define HH_Parameter_Communication_HH

#include "Parameters.hpp"

namespace laplace {

/* Function to let rank 0 communicate the parameters of the problem to the other ranks */
void communicate(parameters& p, int rank, int size);

/* Overload that gets rank and size on its own */
void communicate(parameters& p);

}

#endif