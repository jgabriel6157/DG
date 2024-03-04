#include "Mesh.hxx"

//Constructor definition
Mesh::Mesh(int numCells, double length) : numCells(numCells)
{
    //Generate cells
    for (int i = 0; i < numCells; i++) {
        Cell cell;
        //Calculate length of each cell
        cell.cellLength = length/numCells; //Assuming cells are uniform length
        //Initialize vertices of the cell
        cell.vertices.push_back(i*cell.cellLength); //1D mesh
        //Initialize neighbors of the cell
        cell.neighbors.push_back((i - 1 + numCells) % numCells); //Periodic BC
        cell.neighbors.push_back((i + 1) % numCells); //Periodic BC

        cells.push_back(cell);
    }
}

//Accessor function definition for cells
const std::vector<Cell>& Mesh::getCells() const 
{
    return cells;
}

//Return the number of cells
const int& Mesh::getNumCells() const
{
    return numCells;
}
