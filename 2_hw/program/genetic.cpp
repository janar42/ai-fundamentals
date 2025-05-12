//---------------------------------------------------------------------------
//
//  Template for solving the Travelling Salesman Problem
//  (c) 2021 Ladislava Smítková Janků <ladislava.smitkova@fit.cvut.cz>
//
//  genetic.cpp: Implementation of genetic algorithm with these parameters:
//
//  Population:        500 individuals (can be modified by POPULATION constant)
//  Generations:       30 (can be modified by GENERATIONS constant)
//  Crossover method:  OX or PMX
//  Mutation method:   reversion of the path
//
//  Crossover probability:    95%  (PROBABILITY_CROSSOVER)
//  Mutation probability:     stepped by 5%  (PROBABILITY_MUT_STEP)
//
//  If the fitness value of the actual generation is better than last one,
//  mutation probability will be set to zero. In other case, mutation
//  probability will be increased by PROBABILITY_MUT_STEP.
//
//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <string.h>
#include <limits.h>
#include <algorithm>
#include <limits>
#include <random>
#include <set>
#include <map>
#include <iomanip>

#include "genetic.h"
#include "path.h"

#define POPULATION              100
#define PROBABILITY_CROSSOVER   0.95
#define PROBABILITY_MUT_STEP    0.05
#define GENERATIONS             100

using namespace std;

struct TIndividual {
    std::vector<int> path;
    double fitness;
    double Pr;
    double q;

    bool operator<(const TIndividual & other) const {
        return this->fitness > other.fitness;
    }
};

std::vector<TIndividual> individuals;

bool compareByFitness(const TIndividual &a, const TIndividual &b) {
    return a.fitness > b.fitness;
}

double computeFitness(const std::vector<int> &path, TMatrix *matrix) {
    double fitnessReversed = 0;
    for (int i = 0; i < path.size() - 1; i++) {
        fitnessReversed += matrix->getDistance(path[i], path[i + 1]);
    }
    fitnessReversed += matrix->getDistance(path.front(), path.back());

    return 1 / fitnessReversed;
}

void recalculate(TMatrix *matrix) {
    // calculate fitness and Pr
    double sum_fitness = 0;
    for (auto ind = individuals.begin(); ind != individuals.end(); ++ind) {
        ind->fitness = computeFitness(ind->path, matrix);
        sum_fitness += ind->fitness;
    }

    // sort population by fitness
    std::sort(individuals.begin(), individuals.end(), compareByFitness);

    // compute Pr_i and q_i
    double q = 0;
    for (auto ind = individuals.begin(); ind != individuals.end(); ++ind) {
        ind->Pr = ind->fitness / sum_fitness;
        q += ind->Pr;
        ind->q = q;
    }
}

void selection() {
    vector<TIndividual> newGeneration;

    random_device rd;
    mt19937 gen(rd());

    for (auto i = 0; i < POPULATION; i++) {
        std::set<TIndividual> tournament;

        std::uniform_int_distribution<int> distribution(0, POPULATION - 1);
        for (int j = 0; j < 5; j++) {
            tournament.insert(individuals[distribution(gen)]);
        }
        newGeneration.push_back(*tournament.begin());
    }

    individuals = newGeneration;
}

bool pathContainsCity(const std::vector<int> &path, int city) {
    for (auto x = path.begin(); x != path.end(); ++x) {
        if (city == *x) return true;
    }
    return false;
}

void doCrossoverOX(std::vector<TIndividual> &result, TMatrix *matrix, TIndividual &a, TIndividual &b) {
    TIndividual aa, bb;
    aa = a;
    bb = b;

    size_t size = a.path.size();
    size_t copyFrom = size / 3;
    size_t copyTo = copyFrom * 2;

    set<int> aCenter;
    set<int> bCenter;

    for (size_t i = copyFrom; i < copyTo; i++) {
        aCenter.insert(a.path[i]);
        bCenter.insert(b.path[i]);
    }

    for (size_t i = copyTo + 1, j = copyTo + 1, times = 0; times < size; times++) {
        if (i >= copyFrom && i <= copyTo) {
            i = (i + 1) % size;
            continue;
        }
        if (aCenter.find(b.path[j]) == aCenter.end()) {
            aa.path[i] = b.path[j];
            i = (i + 1) % size;
            j = (j + 1) % size;
        } else {
            j = (j + 1) % size;
        }
    }

    for (size_t i = copyTo + 1, j = copyTo + 1, times = 0; times < size; times++) {
        if (i >= copyFrom && i <= copyTo) {
            i = (i + 1) % size;
            continue;
        }
        if (bCenter.find(a.path[j]) == bCenter.end()) {
            bb.path[i] = a.path[j];
            i = (i + 1) % size;
            j = (j + 1) % size;
        } else {
            j = (j + 1) % size;
        }
    }

    result.push_back(aa);
    result.push_back(bb);
}

void doCrossoverPMX(std::vector<TIndividual> &result, TMatrix *matrix, TIndividual &a, TIndividual &b) {
    TIndividual aa, bb;
    aa = a;
    bb = b;

    size_t size = a.path.size();
    size_t third = size / 3;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distribution(1, size - third - 1);

    size_t copyFrom = distribution(gen);
    size_t copyTo = copyFrom + third;

    set<int> aCenter;
    set<int> bCenter;
    map<int, int> aMap;
    map<int, int> bMap;

    for (size_t i = copyFrom; i <= copyTo; i++) {
        aCenter.insert(a.path[i]);
        bCenter.insert(b.path[i]);
        aMap[a.path[i]] = b.path[i];
        bMap[b.path[i]] = a.path[i];
    }

    for (size_t i = 0; i < size; i++) {
        if (i >= copyFrom && i <= copyTo) continue;

        if (aCenter.find(b.path[i]) == aCenter.end()) {
            aa.path[i] = b.path[i];
        } else {
            int tmp = aMap[b.path[i]];
            while (aCenter.find(tmp) != aCenter.end()) {
                tmp = aMap[tmp];
            }
            aa.path[i] = tmp;
        }
    }

    for (size_t i = 0; i < size; i++) {
        if (i >= copyFrom && i <= copyTo) continue;

        if (bCenter.find(a.path[i]) == bCenter.end()) {
            bb.path[i] = a.path[i];
        } else {
            int tmp = bMap[a.path[i]];
            while (bCenter.find(tmp) != bCenter.end()) {
                tmp = bMap[tmp];
            }
            bb.path[i] = tmp;
        }
    }

    result.push_back(aa);
    result.push_back(bb);
}

void crossover(TMatrix *matrix, TCrossoverMethod crossoverMethod) {
    std::vector<TIndividual> crossedOver;
    std::vector<TIndividual>::iterator candidate = individuals.end();

    for (auto ind = individuals.begin(); ind != individuals.end(); ++ind) {
        // select candidates to the crossover process
        if (drand48() <= PROBABILITY_CROSSOVER) {
            if (candidate == individuals.end()) {
                // this is the first parent
                candidate = ind;
            } else {
                // now we have both parents, we can do a crossover
                if (crossoverMethod == CROSSOVER_METHOD_PMX)
                    doCrossoverPMX(crossedOver, matrix, *ind, *candidate);
                else
                    doCrossoverOX(crossedOver, matrix, *ind, *candidate);

                candidate = individuals.end();
            }
        } else
            crossedOver.push_back(*ind);
    }

    // If we got odd number of parents, do nothing with this candidate and push it directly
    // into the new generation.
    if (candidate != individuals.end())
        crossedOver.push_back(*candidate);

    // crossover is done
    individuals = crossedOver;
}

void mutation(double probability) {
    for (auto indivi = individuals.begin(); indivi != individuals.end(); ++indivi) {
        if (drand48() <= probability) {
            random_device rd;
            mt19937 gen(rd());

            uniform_int_distribution<> distribution(1, indivi->path.size() - 5);
            int from = distribution(gen);

            shuffle(indivi->path.begin() + from, indivi->path.begin() + from + 4, gen);
        }
    }
}

void printState(int generation) {
    cout << "[" << generation << "]  " << std::setprecision(20) << individuals.at(0).fitness * 10000;
    cout << "  {";
    for (int i : individuals.at(0).path) {
        cout << i << " ";
    }
    cout << "}" << endl;
}

std::vector<int> salesmanProblemGenetic(TMatrix *matrix, TCrossoverMethod crossoverMethod) {
    unsigned i, j;
    int x;
    std::vector<int>::iterator p;
    double mutation_probability = 0;
    double lastFitness = -std::numeric_limits<double>::max();
    std::vector<int> best;
    double bestFitness = -std::numeric_limits<double>::max();

    // initialization of random number generator
    srand(getpid());

    // born first population
    for (i = 0; i < POPULATION; i++) {
        TIndividual ind;

        // Generate some random path: Place city indexes to the vector in some random order.
        // At index 0 will be city we start from.
        ind.path.clear();
        ind.path.push_back(0);
        j = 1;
        while (j < matrix->getNumberOfTargets()) {
            x = random() % matrix->getNumberOfTargets();
            p = find(ind.path.begin(), ind.path.end(), x);
            if (p == ind.path.end()) {
                ind.path.push_back(x);
                j++;
            }
        }

        // Store this path into table of individuals.
        // Fitness and other parameters will be computaed later.
        individuals.push_back(ind);
    }

    // compute fitnesses and sort individuals
    recalculate(matrix);
    printState(0);

    // remember the best solution
    best = individuals.at(0).path;
    bestFitness = individuals.at(0).fitness;

    // run simulation
    for (i = 1; i < GENERATIONS; i++) {
        // selection: select individuals for a new generation
        selection();

        // crossover
        crossover(matrix, crossoverMethod);

        // mutation
        if (mutation_probability > 0) mutation(mutation_probability);

        // print the best result
        recalculate(matrix);
        printState(i);

        // if fitness < lastFitness, increase mutation probability by one step
        if (individuals.at(0).fitness < lastFitness)
            mutation_probability += PROBABILITY_MUT_STEP;
        else
            mutation_probability = 0;

        lastFitness = individuals.at(0).fitness;

        if (lastFitness > bestFitness) {
            best = individuals.at(0).path;
            bestFitness = individuals.at(0).fitness;
        }
    }

    return best;
}
