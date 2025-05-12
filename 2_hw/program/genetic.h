//-------------------------------------------------------------------
//
//  Template for solving the Travelling Salesman Problem
//  (c) 2021 Ladislava Smítková Janků <ladislava.smitkova@fit.cvut.cz>
//
//-------------------------------------------------------------------

#ifndef GENETIC_H
#define GENETIC_H

#include "tmatrix.h"
#include <vector>

typedef enum { CROSSOVER_METHOD_OX, CROSSOVER_METHOD_PMX } TCrossoverMethod;

std::vector<int> salesmanProblemGenetic( TMatrix *matrix, TCrossoverMethod crossoverMethod);

#endif // GENETIC_H
