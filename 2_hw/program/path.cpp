//-------------------------------------------------------------------
//
//  Template for solving the Travelling Salesman Problem
//  (c) 2021 Ladislava Smítková Janků <ladislava.smitkova@fit.cvut.cz>
//
//  path.cpp: Misc utils for the work with paths
//
//-------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "path.h"

int calculatePathLength( const std::vector<int> &test, TMatrix *matrix)
{
    int len = 0;
    int last = -1;
    for (auto i=test.begin(); i != test.end(); ++i) {
        if (last != -1) {
            len += matrix->getDistance( last, *i);
        }
        last = *i;
    }
    return len + matrix->getDistance( last, test[0]);
}

void printPath( const std::vector<int> &path)
{
    printf("[");
    for (auto i = path.begin(); i != path.end(); ++i) {
        if (i != path.begin()) printf(" ");
        printf("%d", *i);
    }
    printf(" %d]", path[0]);
}

void printPathWithNames( const std::vector<int> &path, TMatrix *matrix)
{
    bool first = true;
    printf("[");
    for (auto i=path.begin(); i != path.end(); ++i) {
        if (!first) printf(" - ");
        printf("%s", matrix->getTargetName(*i));
        first = false;
    }
    printf(" - %s]", matrix->getTargetName( path[0]));
}
