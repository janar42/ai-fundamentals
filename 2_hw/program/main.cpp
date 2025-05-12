//-------------------------------------------------------------------
//
//  Template for solving the Travelling Salesman Problem
//  (c) 2021 Ladislava Smítková Janků <ladislava.smitkova@fit.cvut.cz>
//
//  main.cpp: This source contains global function main,
//  which is the designated start of the program.
//
//-------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tmatrix.h"
#include "path.h"
#include "bruteforce.h"
#include "genetic.h"

#define ALGO_UNKNOWN      0
#define ALGO_BRUTEFORCE   1
#define ALGO_GENETIC_OX   2
#define ALGO_GENETIC_PMX  3

TMatrix *matrix = NULL;

void usage(char *appname) {
    char *app = strrchr(appname, '/');
    if (app == NULL) app = appname;
    else app++;

    printf("Template for solving the Travelling Salesman Problem.\n"
           "(c) 2021 Ladislava Smítková Janků <ladislava.smitkova@fit.cvut.cz>\n\n"
           "Usage: %s <distances-file-csv> <algorithm>\n\n"
           "<distances-file-csv>   CSV file with distances. Each line contains\n"
           "                       the name of the city and distances to other\n"
           "                       cities separated by commas.\n\n"
           "                       Example of the content of the input CSV file:\n"
           "                           Prague,0,2215,2292\n"
           "                           Madrid,2215,0,3982\n"
           "                           Ankara,2292,3982,0\n\n"
           "                       This means that Madrid is from Ankara 3982 km far,\n"
           "                       Prague from Madrid 2215 km far and Ankara from\n"
           "                       Prague 2292 km far.\n\n"
           "<algorithm>            Algorithm to be used:\n\n"
           "                       bruteforce  Use brute force; try all variants.\n"
           "                       genetic-ox  Use genetic algorithm (OX crossover).\n"
           "                       genetic-pmx Use genetic algorithm (PMX crossover).\n\n", app);
    exit(1);
}

int main(int argc, char *argv[]) {
    std::vector<int> result;

    // Parse arguments.

    if (argc < 3) usage(argv[0]);

    int algorithm = ALGO_UNKNOWN;
    if (strcmp(argv[2], "bruteforce") == 0) algorithm = ALGO_BRUTEFORCE;
    else if (strcmp(argv[2], "genetic-ox") == 0) algorithm = ALGO_GENETIC_OX;
    else if (strcmp(argv[2], "genetic-pmx") == 0) algorithm = ALGO_GENETIC_PMX;
    else
        usage(argv[0]);

    try {
        // Load matrix with distances and run the algorithm.
        TMatrix *matrix = new TMatrix(argv[1]);
        try {
            matrix->printMatrix();

            printf("\nSolving...\n");

            switch (algorithm) {
                case ALGO_BRUTEFORCE:
                    result = salesmanProblemBruteForce(matrix);
                    break;

                case ALGO_GENETIC_OX:
                    result = salesmanProblemGenetic(matrix, CROSSOVER_METHOD_OX);
                    break;

                case ALGO_GENETIC_PMX:
                    result = salesmanProblemGenetic(matrix, CROSSOVER_METHOD_PMX);
                    break;

                default:
                    throw "unknown algorithm";
            }

            printf("\nResult:\n");
            printPathWithNames(result, matrix);
            printf("\nPath length: %d\n", calculatePathLength(result, matrix));
        }
        catch (...) {
            delete matrix;
            throw;
        }
        delete matrix;
    }
    catch (const char *err) {
        fprintf(stderr, "Exception: %s\n", err);
        exit(1);
    }
}
