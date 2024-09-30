#ifndef MESHHEADERDEF
#define MESHHEADERDEF

#include <vector>

struct Cell
{
    std::vector<double> vertices; //Vertex coordinates
    std::vector<int> neighbors; //Indices of neighboring cells
    double dx; //Length of cell in x
};

class Mesh
{
private:
    std::vector<Cell> cells;
    int nx;
    int nvx;
    double domainLengthX;
    double domainMaxVX;
    double dvx; //Distance between velocity points

public:
    //Constructor
    Mesh(int nx, int nvx, double domainLengthX, double domainMaxVX);

    //Accessor
    const std::vector<Cell>& getCells() const;

    const int& getNX() const;

    const int& getNVX() const;

    double getDVX() const;

    double getVelocity(int velocityIndex) const;
};



#endif