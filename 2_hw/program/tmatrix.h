//-------------------------------------------------------------------
//
//  Template for solving the Travelling Salesman Problem
//  (c) 2021 Ladislava Smítková Janků <ladislava.smitkova@fit.cvut.cz>
//
//-------------------------------------------------------------------

#ifndef TMATRIX_H
#define TMATRIX_H

class TMatrix
{
public:
    TMatrix( const char *csvFileName);
    TMatrix( int numberOfTargets_);
    ~TMatrix();

    unsigned getNumberOfTargets();
    const char *getTargetName( int target);
    void setTargetName( int target, const char *name);
    int getDistance( int fromTarget, int toTarget);
    void setDistance( int fromTarget, int toTarget, int distance);
    void printMatrix();

private:
    int numberOfTargets;
    int (*matrix)[];
    char *(*targetName)[];
};

#endif // TMATRIX_H
