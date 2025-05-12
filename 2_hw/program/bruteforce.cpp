//-------------------------------------------------------------------
//
//  Template for solving the Travelling Salesman Problem
//  (c) 2021 Ladislava Smítková Janků <ladislava.smitkova@fit.cvut.cz>
//
//  bruteforce.cpp: Implementation of the Brute Force Algorithm
//
//-------------------------------------------------------------------

#include "bruteforce.h"
#include "limits.h"
#include "path.h"
#include <bits/stdc++.h>

void solveBF( unsigned from, std::vector<int> &test, std::vector<int> &best, int &bestLength, TMatrix *matrix)
{
    if (from == test.size() - 1) {

        // Cannot modify next parts of the vector, we are at the end. So test this vector.

        printf("\rtesting: ");
        printPath( test);
        printf(" ");
        fflush( stdout);

        int len = calculatePathLength( test, matrix);
        if (len < bestLength) {

            // We have new path with better distance. Save it into variable 'best'.
            // Save the path length into the variable 'bestLength'.

            bestLength = len;
            best = test;

            printf("\rnew solution (len=%d): ", bestLength);
            printPath( best);
            printf("\n");
        }
    }
    else {
        unsigned i;
        int tmp;

        // For each city (set into the vector to the index given by variable 'from'),
        // call this procedure recursively with variable 'from' incremented by 1.
        // For these purposes, we must always exchange 2 cities, because use of the city
        // must be unique.

        for (i = from; i < test.size(); i++) {
            // Exchange the city given by index 'from' and the other city in the last part
            // of the vector.
            tmp = test[i];
            test[i] = test[from];
            test[from] = tmp;

            // Try this vector. We can assume that in this recursive call will be tested
            // all combination of all cities contained in the vector from the index 'from+1'
            // to the end of this vector.
            solveBF( from + 1, test, best, bestLength, matrix);

            // Return both cities back to their original positions.
            tmp = test[i];
            test[i] = test[from];
            test[from] = tmp;
        }
    }
}

std::vector<int> salesmanProblemBruteForce( TMatrix *matrix)
{
    std::vector<int> test, best;
    int bestLength = INT_MAX;
    unsigned i;

    // Generate initial sequence of cities.
    // Place city indexes to the vector in some random order.
    // At index 0 will be city we start from.

    for (i=0; i<matrix->getNumberOfTargets(); i++)
        test.push_back(i);

    // Now we can search the best path from city with index 0 through
    // all other cities and back to the starting city.
    // Therefore the first index will be always the same (index of starting city).
    // The index of the last city before we return back to the city we start from is
    // at the end of the vector.
    // Use recursive algorithm - on the top of the stack we will try to put (one by one)
    // all cities to the index 1. Other positions on the vector will be modified later
    // on recursive calls.

    bestLength = INT_MAX;
    solveBF( 1, test, best, bestLength, matrix);
    return best;
}
