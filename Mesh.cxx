#include "Mesh.hxx"
#include <cmath>

//Constructor definition
Mesh::Mesh(int nx, int nvx, double domainLengthX, double domainMaxVX) : nx(nx), nvx(nvx)
{
    //Generate cells
    for (int i = 0; i < nx; i++) 
    {
        for (int j = 0; j < nvx; j++)
        {
            Cell cell;
            //Calculate length of each cell
            cell.dx = domainLengthX/nx; //Assuming cells are uniform length in x
            cell.dvx = 2.0*domainMaxVX/nvx; //Assuming cells are uniform length in vx and -domainMaxVX<vx<domainMaxVX
            //Initialize vertices of the cell
            cell.vertices.push_back({i*cell.dx, -domainMaxVX+j*cell.dvx}); 
            cell.vertices.push_back({(i+1)*cell.dx, -domainMaxVX+j*cell.dvx});
            cell.vertices.push_back({i*cell.dx, -domainMaxVX+(j+1)*cell.dvx});
            cell.vertices.push_back({(i+1)*cell.dx, -domainMaxVX+(j+1)*cell.dvx});
            //Initialize neighbors of the cell
            cell.neighbors.push_back((i-1+nx) % nx); //Periodic BC in x
            cell.neighbors.push_back((i+1) % nx); //Periodic BC in x
            if (j==0) //0 boundary condition in dvx, neighbors should not be called for j=0, j=nvx-1
            {
                cell.neighbors.push_back(0);
                cell.neighbors.push_back(j+1);
            }
            else if (j==nvx-1)
            {
                cell.neighbors.push_back(j-1);
                cell.neighbors.push_back(nvx-1);
            }
            else
            {
                cell.neighbors.push_back(j-1);
                cell.neighbors.push_back(j+1);
            }

            cells.push_back(cell);
        }
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

//Return the number of cells in vx
const int& Mesh::getNVX() const
{
    return nvx;
}
