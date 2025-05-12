//-------------------------------------------------------------------
//
//  Template for solving the Travelling Salesman Problem
//  (c) 2021 Ladislava Smítková Janků <ladislava.smitkova@fit.cvut.cz>
//
//  tmatrix.cpp: Implementation of the TMatrix class.
//
//-------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tmatrix.h"

// Calculate the length of the string with respecting of the UTF-8 encoding.
int strlen_utf8(const char *s) {
    int count = 0;
    while (*s) {
        count += (*s++ & 0xC0) != 0x80;
    }
    return count;
}

// Read the city name from input file.
char *readCityName(char **p) {
    char *from = *p;
    while ((**p != 0) && (**p != ',') && (**p != '\r') && (**p != '\n')) (*p)++;
    if (!**p) throw "CSV syntax error (readCityName)";
    *((*p)++) = 0;
    return from;
}

// Read the distance from the input file.
int readNumber(char **p) {
    int result = 0;

    if (!**p) return -1;

    while ((**p >= '0') && (**p <= '9')) {
        result = result * 10 + (**p) - '0';
        (*p)++;
    }
    if (**p != ',') {
        if (**p == '\r') (*p)++;
        if (**p != '\n') throw "CSV syntax error (readNumber)";
    }
    (*p)++;
    return result;
}

TMatrix::TMatrix(const char *csvFileName) {
    FILE *f;
    char *line = NULL;
    size_t len = sizeof(line);
    char *p;
    int i;
    int cityNum = 0;

    // Open input CSV file and load distance matrix and names of all cities.
    if ((f = fopen(csvFileName, "r")) == NULL)
        throw "Cannot open CSV file.";

    try {
        // Read all lines from the input CSV file.
        while (getline(&line, &len, f) != -1) {

            // Pointer 'p' contains position in line for reading the next field.
            p = line;

            // The first field must be the name of the city.
            // Pointer 'p' will be set to the next field.
            readCityName(&p);

            // Are we here for the first time?
            if (!matrix) {
                // Yes. We need to calculate the number of cities.
                // This value is equal to the number of fields followed after the city name.

                // Backup the value of the pointer 'p'
                char *firstDist = p;

                // Now read the number of fields after the city name. This value is equal
                // to the number of cities.
                numberOfTargets = 0;
                while ((readNumber(&p)) >= 0) numberOfTargets++;

                // We need to allocate memory for city names.
                targetName = (char *(*)[]) malloc(numberOfTargets * sizeof(char *));

                // For this moment, set all names to NULL.
                // Yes, we don't need it, but it's a good thing for debugging.
                bzero(targetName, numberOfTargets * sizeof(char *));

                // It's time to allocate memory for the distance matrix.
                // This matrix has a dimension of numberOfTargets * numberOfTargets.
                matrix = (int (*)[]) malloc(numberOfTargets * numberOfTargets * sizeof(int));

                // It's a good idea to fill this memory by zeroes.
                // This is the same case - we don't need it, but it's a good thing for debugging.
                bzero(matrix, numberOfTargets * numberOfTargets * sizeof(int));

                // Restore the value of the pointer 'p'.
                // Pointer 'p' will be set to the next field after the name of the city.
                p = firstDist;
            }

            // Store the city name.
            setTargetName(cityNum, line);

            // Store all distances into the matrix.
            for (i = 0; i < numberOfTargets; i++) {
                setDistance(cityNum, i, readNumber(&p));
            }

            // It's done. Increment the variable 'cityNum', which contains an index of actual
            // loaded city.
            cityNum++;
        }

        // free the memory allocated by function getline().
        free(line);
    }
    catch (...) {
        fclose(f);
        throw;
    }
    fclose(f);
}

TMatrix::TMatrix(int numberOfTargets_) {
    // This constructor can be used when we wish to create the distance matrix manually.
    // In this case, you can fill allocated memory with zeroes by yourself.
    // If you wish it, it's not necessary.
    numberOfTargets = numberOfTargets_;
    targetName = (char *(*)[]) malloc(numberOfTargets_ * sizeof(char *));
    matrix = (int (*)[]) malloc(numberOfTargets_ * numberOfTargets_ * sizeof(int));
}

TMatrix::~TMatrix() {
    // Frees all allocated memory for the matrix and city names.
    if (matrix) {
        free(matrix);
        matrix = NULL;
    }
    if (targetName) {
        int i;
        for (i = 0; i < numberOfTargets; i++) free((*targetName)[i]);
        free(targetName);
        targetName = NULL;
    }
    numberOfTargets = 0;
}

unsigned TMatrix::getNumberOfTargets() {
    return numberOfTargets;
}

const char *TMatrix::getTargetName(int target) {
    if ((target < 0) || (target >= numberOfTargets)) throw "TMatrix::getTargetName: Bad target.";
    return (*targetName)[target];
}

void TMatrix::setTargetName(int target, const char *name) {
    if ((target < 0) || (target >= numberOfTargets)) throw "TMatrix::setTargetName: Bad target.";
    (*targetName)[target] = strdup(name);
}

void TMatrix::setDistance(int fromTarget, int toTarget, int distance) {
    if (distance < 0) throw "TMatrix::setDistance: Negative value of the distance.";
    if ((fromTarget < 0) || (fromTarget >= numberOfTargets)) throw "TMatrix::setDistance: Bad target (fromTarget).";
    if ((toTarget < 0) || (toTarget >= numberOfTargets)) throw "TMatrix::setDistance: Bad target (toTarget).";

    (*matrix)[fromTarget * numberOfTargets + toTarget] = distance;
    (*matrix)[toTarget * numberOfTargets + fromTarget] = distance;
}

int TMatrix::getDistance(int fromTarget, int toTarget) {
    if ((fromTarget < 0) || (fromTarget >= numberOfTargets)) throw "TMatrix::getDistance: Bad target (fromTarget).";
    if ((toTarget < 0) || (toTarget >= numberOfTargets)) throw "TMatrix::getDistance: Bad target (toTarget).";

    return (*matrix)[fromTarget * numberOfTargets + toTarget];
}

void TMatrix::printMatrix() {
    int i, j;
    for (i = 0; i < numberOfTargets; i++) {
        printf("%s%*s", (*targetName)[i], 10 - strlen_utf8((*targetName)[i]), "");
        for (j = 0; j < numberOfTargets; j++)
            printf(" %4d", getDistance(i, j));
        printf("\n");
    }
}
