#include "Mesh.hxx"

//Constructor definition
Mesh::Mesh(int nvx, double domainMaxVX) : nvx(nvx)
{
    //Generate cells
    for (int j = 0; j < nvx; j++) {
        Cell cell;
        //Calculate length of each cell
        cell.dvx = 2.0*domainMaxVX/nvx; //Assuming cells are uniform length
        //Initialize vertices of the cell
        cell.vertices.push_back(-domainMaxVX+j*cell.dvx); //1D mesh
        cell.vertices.push_back(-domainMaxVX+(j+1)*cell.dvx);
        //Initialize neighbors of the cell

        if (j==0) //0 boundary condition in dvx, neighbors should not be called for j=0, j=nvx-1
        {
            cell.neighbors.push_back(0); //Bottom Neighbor
            cell.neighbors.push_back((j+1)); //Top Neighbor
        }
        else if (j==nvx-1)
        {
            cell.neighbors.push_back((j-1)); //Bottom Neighbor
            cell.neighbors.push_back((nvx-1)); //Top Neighbor
        }
        else
        {
            cell.neighbors.push_back((j-1)); //Bottom Neighbor
            cell.neighbors.push_back((j+1)); //Top Neighbor
        }

        cells.push_back(cell);
    }
}

//Accessor function definition for cells
const std::vector<Cell>& Mesh::getCells() const 
{
    return cells;
}

//Return the number of cells
const int& Mesh::getNVX() const
{
    return nvx;
}
