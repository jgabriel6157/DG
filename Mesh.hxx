#ifndef MESHHEADERDEF
#define MESHHEADERDEF

#include <vector>

struct Cell
{
    std::vector<double> vertices; //Vertex coordinates
    std::vector<int> neighbors; //Indices of neighboring cells
    double dvx; //Length of cell
};

class Mesh
{
private:
    std::vector<Cell> cells;
    int nvx;

public:
    //Constructor
    Mesh(int nvx, double domainMaxVX);

    //Accessor
    const std::vector<Cell>& getCells() const;

    const int& getNVX() const;
};



#endif