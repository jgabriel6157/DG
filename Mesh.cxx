#include "Mesh.hxx"
#include <cmath>

/*
                             N3
                             ^
                   v2________|__________v3
                    |        |f3        |
                    |                   |
                    |f0               f1|
              N0 <--|--     cell      --|--> N1
                    |                   |
                    |                   |                    
                    |________|f2________|
                   v0        |          v1 
                             v  
                             N2   
*/


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
            // cell.dvx = 2.0*domainMaxVX/nvx; //Assuming cells are uniform length in vx and -domainMaxVX<vx<domainMaxVX
            // //Initialize vertices of the cell 
            // cell.vertices.push_back({i*cell.dx, -domainMaxVX+j*cell.dvx}); //Bottom Left 
            // cell.vertices.push_back({(i+1)*cell.dx, -domainMaxVX+j*cell.dvx}); //Bottom Right
            // cell.vertices.push_back({i*cell.dx, -domainMaxVX+(j+1)*cell.dvx}); //Top Left
            // cell.vertices.push_back({(i+1)*cell.dx, -domainMaxVX+(j+1)*cell.dvx}); //Top Right
            cell.dvx = 1.0/nvx; //Assuming cells are uniform length in vx and -domainMaxVX<vx<domainMaxVX
            //Initialize vertices of the cell 
            cell.vertices.push_back({i*cell.dx, domainMaxVX+j*cell.dvx}); //Bottom Left 
            cell.vertices.push_back({(i+1)*cell.dx, domainMaxVX+j*cell.dvx}); //Bottom Right
            cell.vertices.push_back({i*cell.dx, domainMaxVX+(j+1)*cell.dvx}); //Top Left
            cell.vertices.push_back({(i+1)*cell.dx, domainMaxVX+(j+1)*cell.dvx}); //Top Right
            //Initialize neighbors of the cell {Left, Right, Bottom, Top}
            cell.neighbors.push_back(((i-1+nx) % nx)+j*nx); //Periodic BC in x (Left Neighbor)
            cell.neighbors.push_back(((i+1) % nx)+j*nx); //Periodic BC in x (Right Neighbor)
            if (j==0) //0 boundary condition in dvx, neighbors should not be called for j=0, j=nvx-1
            {
                cell.neighbors.push_back(i+0*nx); //Bottom Neighbor
                cell.neighbors.push_back(i+(j+1)*nx); //Top Neighbor
            }
            else if (j==nvx-1)
            {
                cell.neighbors.push_back(i+(j-1)); //Bottom Neighbor
                cell.neighbors.push_back(i+(nvx-1)); //Top Neighbor
            }
            else
            {
                cell.neighbors.push_back(i+(j-1)); //Bottom Neighbor
                cell.neighbors.push_back(i+(j+1)); //Top Neighbor
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
