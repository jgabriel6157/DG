#include "Mesh.hxx"

//Constructor definition
Mesh::Mesh(int nx, int nvx, double domainLengthX, double domainMaxVX) : nx(nx), nvx(nvx), domainLengthX(domainLengthX), domainMaxVX(domainMaxVX)
{
    //Generate cells
    for (int i = 0; i < nx; i++) 
    {
        Cell cell;
        //Calculate length of each cell
        cell.dx = domainLengthX/nx; //Assuming cells are uniform length
        //Initialize vertices of the cell
        cell.vertices.push_back(i*cell.dx); //1D mesh
        //Initialize neighbors of the cell
        cell.neighbors.push_back((i - 1 + nx) % nx); //Periodic BC
        cell.neighbors.push_back((i + 1) % nx); //Periodic BC

        cells.push_back(cell);
    }
}

//Accessor function definition for cells
const std::vector<Cell>& Mesh::getCells() const 
{
    return cells;
}

//Return the number of cells in x
const int& Mesh::getNX() const
{
    return nx;
}

//Return the number of cells in x
const int& Mesh::getNVX() const
{
    return nvx;
}

double Mesh::getVelocity(int velocityIndex) const
{
    double dvx = 2.0*domainMaxVX/(nvx-1.0);
    return -domainMaxVX + velocityIndex*dvx;
}
