//-------------------------------------------------------------------
//
//  Template for solving the Travelling Salesman Problem
//  (c) 2021 Ladislava Smítková Janků <ladislava.smitkova@fit.cvut.cz>
//
//-------------------------------------------------------------------

#ifndef PATH_H
#define PATH_H

#include "tmatrix.h"
#include <vector>

int calculatePathLength( const std::vector<int> &test, TMatrix *matrix);
void printPath( const std::vector<int> &path);
void printPathWithNames( const std::vector<int> &path, TMatrix *matrix);

#endif // PATH_H
