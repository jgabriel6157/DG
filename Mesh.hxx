#ifndef MESHHEADERDEF
#define MESHHEADERDEF

#include <vector>

struct Cell
{
    std::vector<double> vertices; //Vertex coordinates
    std::vector<int> neighbors; //Indices of neighboring cells
    double cellLength; //Length of cell
};

class Mesh
{
private:
    std::vector<Cell> cells;
    int numCells;

public:
    //Constructor
    Mesh(int numCells, double length);

    //Accessor
    const std::vector<Cell>& getCells() const;

    const int& getNumCells() const;
};



#endif