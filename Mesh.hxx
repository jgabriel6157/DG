#ifndef MESHHEADERDEF
#define MESHHEADERDEF

#include <vector>

struct Cell
{
    std::vector<double> vertices; //Vertex coordinates
    std::vector<int> neighbors; //Indices of neighboring cells
    double dx; //Length of cell in x
    double dvx; //Distance between velocity points
};

class Mesh
{
private:
    std::vector<Cell> cells;
    int nx;
    int nvx;
    double domainLengthX;
    double domainMaxVX;

public:
    //Constructor
    Mesh(int nx, int nvx, double domainLengthX, double domainMaxVX);

    //Accessor
    const std::vector<Cell>& getCells() const;

    const int& getNX() const;

    const int& getNVX() const;

    double getVelocity(int velocityIndex) const;
};



#endif