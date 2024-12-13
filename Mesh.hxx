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
    int nvy;
    int nvz;
    double domainLengthX;
    double domainMaxVX;
    double domainMaxVY;
    double domainMaxVZ;
    int bc;
    double dvx; //Distance between velocity points
    double dvy;
    double dvz;

public:
    //Constructor
    Mesh(int nx, int nvx, int nvy, int nvz, double domainLengthX, double domainMaxVX, double domainMaxVY, double domainMaxVZ, int bc);

    //Accessor
    const std::vector<Cell>& getCells() const;

    const int& getNX() const;

    const int& getNVX() const;

    const int& getNVY() const;

    const int& getNVZ() const;

    double getDVX() const;

    double getDVY() const;

    double getDVZ() const;

    double getVelocityX(int velocityIndex) const;

    double getVelocityY(int velocityIndex) const;

    double getVelocityZ(int velocityIndex) const;
};



#endif