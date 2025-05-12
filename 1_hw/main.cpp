#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <random>
#include <algorithm>
#include <thread>
#include <stack>

struct Point {
    int x, y, dist;

    Point(int x, int y, int d) : x(x), y(y), dist(d) {};

    Point(int x, int y) : x(x), y(y), dist(0) {};

    Point() : x(0), y(0), dist(0) {};

    bool operator==(const Point &other) const {
        return x == other.x && y == other.y;
    };

    double distance(const Point &target) const {
        return std::sqrt(std::pow((x - target.x), 2) + std::pow((y - target.y), 2));
    };
};

struct heuGreedySearch {
    Point target;

    heuGreedySearch(const Point &target) : target(target) {};

    bool operator()(const Point &p1, const Point &p2) const {
        return p1.distance(target) > p2.distance(target);
    };
};

struct heuAStar {
    Point target;

    heuAStar(const Point &target) : target(target) {};

    bool operator()(const Point &p1, const Point &p2) const {
        return (p1.distance(target) + p1.dist) > (p2.distance(target) + p2.dist);
    };
};

class Maze {
    std::vector<std::vector<char>> maze;
    std::vector<std::vector<Point>> P;
    std::vector<std::vector<int>> D;

    Point start, end;
    size_t expanded = 0;

public:
    Maze(std::ifstream &inputFile) {
        if (!inputFile.is_open()) {
            std::cout << "Soubor nelze otevrit." << std::endl;
            exit(1);
        }

        std::string line, word, comma;

        while (std::getline(inputFile, line)) {
            if (line[0] == 's') {
                std::istringstream in(line);
                in >> word >> start.x >> comma >> start.y;
                in.clear();
            } else if (line[0] == 'e') {
                std::istringstream in(line);
                in >> word >> end.x >> comma >> end.y;
                in.clear();
            } else {
                std::vector<char> row;
                for (char col: line)
                    row.push_back(col);
                maze.push_back(row);
            }
        }
        inputFile.close();
    };

    void Solve(const int alg) {
        if (start == end)
            std::cout << "\nStart a cil jsou na stejnem miste.\n" << std::endl;
        else if (alg == 1)
            BFS();
        else if (alg == 2)
            DFS();
        else if (alg == 3)
            randomSearch();
        else if (alg == 4)
            greedySearch();
        else if (alg == 5)
            AStar();
        else
            std::cout << "Nespravny vstup." << std::endl;
    };

    void printEnd() {
        std::vector<std::vector<char>> print(maze.size(), std::vector<char>(maze[0].size(), ' '));
        for (int i = 0; i < (int) maze.size(); i++) {
            for (int j = 0; j < (int) maze[i].size(); j++) {
                if (maze[i][j] == 'X')
                    print[i][j] = 'X';
                else if (D[i][j] == 0)
                    print[i][j] = 'S';
                else if (end.y == i && end.x == j)
                    print[i][j] = 'E';
                else if (D[i][j] > 0)
                    print[i][j] = '#';
                else
                    print[i][j] = ' ';
            }
        }

        Point curr = {P[end.y][end.x].x, P[end.y][end.x].y};
        while (true) {
            if (curr.x == start.x && curr.y == start.y)
                break;

            print[curr.y][curr.x] = 'o';
            curr = {P[curr.y][curr.x].x, P[curr.y][curr.x].y};
        }

        std::cout << '\n';

        for (size_t i = 0; i < print.size(); i++) {
            for (size_t j = 0; j < print[i].size(); j++) {
                std::cout << print[i][j];
            }
            std::cout << std::endl;
        }

        std::cout
                << "-----------------------------------------\n"
                << "S Start\nE End\n# Opened node\no Path\nX Wall\nspace Fresh node\n"
                << "-----------------------------------------\n"
                << "Nodes expanded: " << expanded << "\n"
                << "Path length: " << D[end.y][end.x]
                << std::endl;
    };

private:
    bool isValid(const int &x, const int &y) {
        return x >= 0 && y >= 0 && maze[y][x] != 'X';
    };

    void BFS() {
        D = std::vector<std::vector<int>>(maze.size(), std::vector<int>(maze[0].size(), -1));
        P = std::vector<std::vector<Point>>(maze.size(), std::vector<Point>(maze[0].size(), {-1, -1}));

        std::queue<Point> Q;
        Q.push(start);
        D[start.y][start.x] = 0;

        int moveX[] = {1, -1, 0, 0};
        int moveY[] = {0, 0, 1, -1};

        while (!Q.empty()) {
            Point curr = Q.front();
            Q.pop();

            for (int i = 0; i < 4; i++) {
                int nextX = curr.x + moveX[i];
                int nextY = curr.y + moveY[i];

                if (isValid(nextX, nextY)) {
                    if (D[nextY][nextX] == -1) {
                        Q.push({nextX, nextY});
                        expanded++;
                        D[nextY][nextX] = D[curr.y][curr.x] + 1;
                        P[nextY][nextX] = curr;

                        if (end.x == nextX && end.y == nextY) {
                            printEnd();
                            return;
                        }
                    }
                }
            }
        }
    };

    void DFS() {
        D = std::vector<std::vector<int>>(maze.size(), std::vector<int>(maze[0].size(), -1));
        P = std::vector<std::vector<Point>>(maze.size(), std::vector<Point>(maze[0].size(), {-1, -1}));

        std::stack<Point> stack;
        stack.push(start);
        D[start.y][start.x] = 0;

        int moveX[] = {1, -1, 0, 0};
        int moveY[] = {0, 0, 1, -1};

        while (!stack.empty()) {
            Point curr = stack.top();
            stack.pop();

            if (curr.x == end.x && curr.y == end.y) {
                printEnd();
                return;
            }

            for (int i = 0; i < 4; i++) {
                int nextX = curr.x + moveX[i];
                int nextY = curr.y + moveY[i];
                if (isValid(nextX, nextY) && D[nextY][nextX] == -1) {
                    stack.push({nextX, nextY});
                    expanded++;
                    D[nextY][nextX] = D[curr.y][curr.x] + 1;
                    P[nextY][nextX] = curr;
                }
            }
        }

    };

    void randomSearch() {
        D = std::vector<std::vector<int>>(maze.size(), std::vector<int>(maze[0].size(), -1));
        P = std::vector<std::vector<Point>>(maze.size(), std::vector<Point>(maze[0].size(), {-1, -1}));

        std::vector<Point> open;
        open.push_back(start);
        D[start.y][start.x] = 0;

        int moveX[] = {1, -1, 0, 0};
        int moveY[] = {0, 0, 1, -1};

        while (!open.empty()) {
            std::random_device rd;
            std::mt19937 rng(rd());
            std::uniform_int_distribution<int> distribution(0, open.size() - 1);
            int rdIndex = distribution(rng);

            Point curr = open[rdIndex];
            open.erase(open.begin() + rdIndex);

            if (curr.x == end.x && curr.y == end.y) {
                printEnd();
                return;
            }
            for (int i = 0; i < 4; i++) {
                int nextX = curr.x + moveX[i];
                int nextY = curr.y + moveY[i];
                Point next = {nextX, nextY};

                if (isValid(nextX, nextY) && D[nextY][nextX] == -1) {
                    open.push_back(next);
                    expanded++;
                    D[nextY][nextX] = D[curr.y][curr.x] + 1;
                    P[nextY][nextX] = curr;
                }
            }
        }
    };

    void greedySearch() {
        D = std::vector<std::vector<int>>(maze.size(), std::vector<int>(maze[0].size(), -1));
        P = std::vector<std::vector<Point>>(maze.size(), std::vector<Point>(maze[0].size(), {-1, -1}));

        std::priority_queue<Point, std::vector<Point>, heuGreedySearch> queue(heuGreedySearch(this->end));
        queue.push(start);
        D[start.y][start.x] = 0;

        int moveX[] = {1, -1, 0, 0};
        int moveY[] = {0, 0, 1, -1};

        while (!queue.empty()) {
            Point curr = queue.top();
            queue.pop();

            for (int i = 0; i < 4; i++) {
                int nextX = curr.x + moveX[i];
                int nextY = curr.y + moveY[i];

                if (end.x == nextX && end.y == nextY) {
                    printEnd();
                    return;
                }

                if (isValid(nextX, nextY) && D[nextY][nextX] == -1) {
                    queue.push({nextX, nextY});
                    expanded++;
                    D[nextY][nextX] = D[curr.y][curr.x] + 1;
                    P[nextY][nextX] = curr;
                }
            }
        }
    };

    void AStar() {
        D = std::vector<std::vector<int>>(maze.size(), std::vector<int>(maze[0].size(), -1));
        P = std::vector<std::vector<Point>>(maze.size(), std::vector<Point>(maze[0].size(), {-1, -1}));

        std::priority_queue<Point, std::vector<Point>, heuAStar> queue(heuAStar(this->end));
        queue.push(start);
        D[start.y][start.x] = 0;

        int moveX[] = {1, -1, 0, 0};
        int moveY[] = {0, 0, 1, -1};

        while (!queue.empty()) {
            Point curr = queue.top();
            queue.pop();

            for (int i = 0; i < 4; i++) {
                int nextX = curr.x + moveX[i];
                int nextY = curr.y + moveY[i];

                if (end.x == nextX && end.y == nextY) {
                    printEnd();
                    return;
                }

                if (isValid(nextX, nextY) && D[nextY][nextX] == -1) {
                    queue.push({nextX, nextY, D[curr.y][curr.x] + 1});
                    expanded++;
                    D[nextY][nextX] = D[curr.y][curr.x] + 1;
                    P[nextY][nextX] = curr;
                }
            }
        }
    };
};

int main(int argc, char *argv[]) {
    int alg;
    std::cout << "Vyberte algoritmus: \n1 - BFS, 2 - DFS, 3 - random search, 4 - greedy search, 5 - A* " << std::endl;
    std::cin >> alg;

    std::ifstream input(argv[1]);

    Maze maze = Maze(input);
    maze.Solve(alg);

    return 0;
}