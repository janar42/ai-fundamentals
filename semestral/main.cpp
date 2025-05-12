#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <filesystem>
#include <random>
#include <set>
#include <map>
#include <algorithm>
#include <limits>
#include <chrono>

#include "src.h"

using namespace std;

#define MAX_ITERS 100000
#define MAX_RESTARTS 100

class Grid {
public:
    size_t size;
    size_t square;
    vector<vector<int>> values;
    vector<vector<bool>> movable;

    Grid() : size(0), square(0), values({}), movable({}) {}

    Grid(const vector<vector<int>> &matrix) {
        size = matrix.size();
        square = sqrt(size);
        movable.resize(size, vector<bool>(size, true));
        values = matrix;

        for (size_t i = 0; i < size; i++) {
            for (size_t j = 0; j < size; j++) {
                if (matrix[i][j] != 0) {
                    movable[i][j] = false;
                }
            }
        }
    }

    void print() const {
        for (size_t i = 0; i < size; i++) {

            for (size_t j = 0; j < size; j++) {
                cout << values[i][j];

                if ((j + 1) % square == 0 && j + 1 != size)
                    cout << "  ";
                else if (j + 1 != size)
                    cout << " ";
            }
            cout << "\n";


            if ((i + 1) % square == 0 && i + 1 != size)
                cout << "\n";
        }
    }

};

class HillClimber {
    int attempts = 0;
    int restarts = 0;
    int minFitness = INT32_MAX;

    Grid starter;
    Grid filled;
    Grid solution;

    bool solved = false;
    set<vector<vector<int>>> visited;

public:

    HillClimber(const Grid & grid) {
        starter = grid;
        filled = {};
        solution = {};
    }

    void fillMatrix() {
        filled = starter;

        for (int i = 0; i < filled.size; i++) {
            set<int> used;

            for (int j = 0; j < filled.size; j++) {
                if (filled.values[i][j] != 0)
                    used.insert(filled.values[i][j]);
            }

            for (int j = 0; j < filled.size; j++) {
                if (filled.values[i][j] == 0) {
                    for (int k = 1; k <= filled.size; k++) {
                        if (!used.count(k)) {
                            filled.values[i][j] = k;
                            used.insert(k);
                            break;
                        }
                    }
                }
            }
        }

        std::cout << "\nDoplnena mrizka:" << std::endl;
        filled.print();
    }

    int antiFitness(Grid &grid) {
        int idealSum = 0;
        for (int i = 1; i <= grid.size; i++) {
            idealSum += i;
        }
        int fitness = 0;

        //counting rows
        for (int i = 0; i < grid.size; i++) {
            int currSum = 0;
            vector<int> count(grid.size + 1, 0);

            for (int j = 0; j < grid.size; j++) {
                currSum += grid.values[i][j];
                count[grid.values[i][j]]++;
            }

            fitness += abs(currSum - idealSum);
            for (int k = 1; k <= grid.size; k++) {
                if (count[k] > 1) { fitness += 1000; }
            }
        }

        //counting columns
        for (int j = 0; j < grid.size; j++) {
            int currSum = 0;
            vector<int> count(grid.size + 1, 0);

            for (int i = 0; i < grid.size; i++) {
                currSum += grid.values[i][j];
                count[grid.values[i][j]]++;
            }

            fitness += abs(currSum - idealSum);
            for (int k = 1; k <= grid.size; k++) {
                if (count[k] > 1) { fitness += 1000; }
            }
        }

        //counting subgrids
        for (int bi = 0; bi < grid.square; bi++) {
            for (int bj = 0; bj < grid.square; bj++) {
                int currSum = 0;
                vector<int> count(grid.size + 1, 0);

                for (int i = 0; i < grid.square; i++) {
                    for (int j = 0; j < grid.square; j++) {
                        int row = bi * grid.square + i;
                        int col = bj * grid.square + j;
                        currSum += grid.values[row][col];
                        count[grid.values[row][col]]++;
                    }
                }

                fitness += abs(currSum - idealSum);
                for (int k = 1; k <= grid.size; k++) {
                    if (count[k] > 1) { fitness += 1000; }
                }
            }
        }

        return fitness;
    }

    bool isCorrectSudoku(Grid &grid) {
        //checking rows
        for (int i = 0; i < grid.size; i++) {
            set<int> used;
            for (int j = 0; j < grid.size; j++) {
                used.insert(grid.values[i][j]);
            }
            if (used.size() != grid.size)
                return false;
        }

        //checking columns
        for (int j = 0; j < grid.size; j++) {
            set<int> used;
            for (int i = 0; i < grid.size; i++) {
                used.insert(grid.values[i][j]);
            }
            if (used.size() != grid.size)
                return false;
        }

        //checking subgrids
        for (int bi = 0; bi < grid.square; bi++) {
            for (int bj = 0; bj < grid.square; bj++) {
                set<int> used;
                for (int i = 0; i < grid.square; i++) {
                    for (int j = 0; j < grid.square; j++) {
                        used.insert(grid.values[(bi * grid.square) + i][(bj * grid.square) + j]);
                    }
                }
                if (used.size() != grid.size)
                    return false;
            }
        }

        return true;
    }

    void swapCells(Grid &grid) {
        static mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
        uniform_int_distribution<int> dist(0, grid.size - 1);
        int a = 0;
        int b = 0;
        int row = 0;

        while (true) {
            a = dist(rng);
            b = dist(rng);
            row = dist(rng);

            if (a == b || !grid.movable[row][a] || !grid.movable[row][b]) {
                continue;
            }

            std::swap(grid.values[row][a], grid.values[row][b]);
            break;
        }
    }

    void shuffleRow(Grid & grid) {
        static mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
        uniform_int_distribution<int> rowDist(0, grid.size - 1);
        uniform_int_distribution<int> numDist(1, grid.size);

        set<int> used;
        int row = rowDist(rng);

        for (int i = 0; i < grid.size; i++) {
            if (starter.values[row][i] != 0)
                used.insert(starter.values[row][i]);
        }

        for (int i = 0; i < grid.size; i++) {
            if (!grid.movable[row][i]) { continue; }

            while (true) {
                int randomNum = numDist(rng);
                if (!used.count(randomNum)) {
                    used.insert(randomNum);
                    grid.values[row][i] = randomNum;
                    break;
                }
            }
        }
    }

    Grid newNeighbor(Grid & grid) {
        mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
        uniform_int_distribution<int> dist01(0, 1);

        Grid currSol = grid;
        swapCells(currSol);

        int random = dist01(rng);
        if (random) {
            shuffleRow(currSol);
        }
        return currSol;
    }

    void execute() {
        fillMatrix();

        while (restarts < MAX_RESTARTS) {
            Grid currSol = filled;
            int currFitness = antiFitness(currSol);
            if (currFitness == 0) {
                solution = currSol;
                solved = true;
                return;
            }
            visited.insert(currSol.values);
            int currTries = 0;

            while (currTries < MAX_ITERS) {
                Grid randNeighbour = newNeighbor(currSol);

                if (visited.count(randNeighbour.values) == 0) {
                    currTries++;

                    currFitness = antiFitness(randNeighbour);
                    currSol = randNeighbour;
                    visited.insert(randNeighbour.values);

                    if (currFitness < minFitness && currFitness != 0) {
                        minFitness = currFitness;
                        solution = randNeighbour;
                        cout << "\n\nNove fitness minimum: " << minFitness << ", Poradi pokusu: " << attempts + currTries << endl;
                        solution.print();
                    }

                    if (currFitness == 0 && isCorrectSudoku(randNeighbour)) {
                        solution = randNeighbour;
                        solved = true;
                        attempts += currTries;
                        return;
                    }
                }
            }
            restarts++;
            if (restarts != MAX_RESTARTS) {
                attempts += currTries;
            }
        }
    }

    void endMessage() {
        if (solved) {
            cout << "\nReseni bylo nalezeno po " << attempts << " pokusech." << endl;
            cout << "Tady je:" << endl;
            solution.print();
        } else {
            cout << "\n\nReseni nebylo nalezeno ani po " << attempts << " pokusech." << endl;
        }
    }

};

int main() {
    int dim, diff;
    cout << "Zadejte rozmer sudoku (4, 9, 16)" << endl;
    cin >> dim;
    if (dim != 4) {
        cout << "Zadejte obtiznost\n1 - light, 2 - medium, 3 - hard" << endl;
        cin >> diff;
    }

    Grid grid;

    if (dim == 9) {
        switch (diff) {
            case 1:
                grid = Grid(in9light);
                break;
            case 2:
                grid = Grid(in9medium);
                break;
            case 3:
                grid = Grid(in9hard);
                break;
            default:
                std::cout << "Spatny vstup" << std::endl;
                return 0;
        }

    } else if (dim == 16) {
        switch (diff) {
            case 1:
                grid = Grid(in16light);
                break;
            case 2:
                grid = Grid(in16medium);
                break;
            case 3:
                grid = Grid(in16hard);
                break;
            default:
                std::cout << "Spatny vstup" << std::endl;
                return 0;
        }
    } else if (dim == 4) {
        grid = Grid(in4);
    } else {
        cout << "Spatny vstup" << endl;
    }

    cout << "\nVybrali jste si tuto sudoku:" << endl;
    grid.print();
    cout << "\n";

    HillClimber hc(grid);
    hc.execute();
    hc.endMessage();


    return 0;
}
