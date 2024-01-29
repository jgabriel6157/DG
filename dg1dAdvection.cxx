#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>

double LegendreP(int n, double x);
double LegendreP_derivative(int n, double x);
double LegendrePorthonormal(int n, double x);
double LegendrePorthonormal_derivative(int n, double x);

void LeastSquares(double uInitialize[], double xj, double dx, int order);
void MatrixMultiply(double M1[20][10], double M2[10][20], double M[20][20], int order);
void MatrixInverse(double M[20][20], double M_inverse[20][20], int order);
void MatrixMultiply2(double M1[20][20], double M2[20][10], double M[20][10], int order);
void MatrixMultiply3(double M1[20][20], double M2[20][20], double M[20][20], int order);

void roots_initial(int n, double x[]);
void roots_final(int n, double x[]);
void weights(int n, double x[], double w[]);
double integrate(int n_func1, int n_func2, int n, double x[], double w[]);
double integrateS(int n_func1, int n_func2, int n, double x[], double w[]);


int main(int argc, char* argv[])
{
    double const pi = M_PI;
    double a = 2*pi;
    int jMax = 5;
    int lMax = 5;
    int tMax = 100000;
    int quadratureOrder = 20;
    double dx = 2*pi/jMax;
    double dt = 0.00001;
    double M[20][20], M_inv[20][20], S[20][20], F1[20][20], F2[20][20];
    // double M_inv[20][20] = {{1.0/dx,0},{0,3.0/dx}};
    // double S[20][20] = {{0,0},{2.0,0}};
    // double F1[20][20] = {{1.0,1.0},{1.0,1.0}};
    // double F2[20][20] = {{1.0,1.0},{-1.0,-1.0}};
    double M_invS[20][20],M_invF1[20][20],M_invF2[20][20];
    double uPre[lMax][jMax];
    double uPost[lMax][jMax];
    double uInitialize[lMax];
    double u1[lMax], u2[lMax];
    double x_roots[quadratureOrder], w[quadratureOrder];
    double xj, val;
    double F_minus,F_plus;
    std::ofstream write_output("Output.csv");
    assert(write_output.is_open());

    roots_initial(quadratureOrder,x_roots);
    roots_final(quadratureOrder,x_roots);
    weights(quadratureOrder,x_roots,w);

    for (int i=0; i<lMax; i++)
    {
        for (int j=0; j<lMax; j++)
        {
            M[i][j] = integrate(i,j,quadratureOrder,x_roots,w)/2.0;
            S[i][j] = integrateS(i,j,quadratureOrder,x_roots,w);
            F1[i][j] = LegendreP(i,1)*LegendreP(j,1);
            F2[i][j] = LegendreP(i,-1)*LegendreP(j,1);
            // if (M[i][j] < 1e-10)
            // {
            //     M[i][j] = 0;
            // }
        }
    }
    
    MatrixInverse(M,M_inv,lMax);

    MatrixMultiply3(M_inv,S,M_invS,lMax);
    MatrixMultiply3(M_inv,F1,M_invF1,lMax);
    MatrixMultiply3(M_inv,F2,M_invF2,lMax);

    for (int j=0; j<jMax; j++)
    {
        xj = j*dx+dx/2.0;
        for (int l=0; l<lMax; l++)
        {
            uInitialize[l] = 0;
        }
        LeastSquares(uInitialize, xj, dx, lMax);
        for (int l=0; l<lMax; l++)
        {
            val = uInitialize[l];
            uPre[l][j] = val;
        }
    }
    for (int t=0; t<tMax; t++)
    {
        for (int j=0; j<jMax; j++)
        {
            for (int l=0; l<lMax; l++)
            {
                uPost[l][j]=0;
                for (int i=0; i<lMax; i++)
                {
                    uPost[l][j]+=M_invS[l][i]*uPre[i][j];
                    uPost[l][j]-=M_invF1[l][i]*uPre[i][j];
                    uPost[l][j]+=M_invF2[l][i]*uPre[i][j-1];
                }
                uPost[l][j]*=a;
                uPost[l][j]*=dt;
                uPost[l][j]/=dx;
                uPost[l][j]+=uPre[l][j];
            }
            // uPost[0][j]=0;
            // uPost[1][j]=0;
            // F_plus = uPre[0][j]+uPre[1][j];
            // F_minus = uPre[0][j-1]+uPre[1][j-1];
            // uPost[0][j] = uPre[0][j]-(a*dt/dx)*(F_plus-F_minus);
            // uPost[1][j] = uPre[1][j]-3.0*(a*dt/dx)*(F_plus+F_minus)+6.0*(a*dt/dx)*uPre[0][j];
        }
    
        for (int j=0; j<jMax; j++)
        {
            for (int l=0; l<lMax; l++)
            {
                uPre[l][j] = uPost[l][j];
            }
        }
    }

    for (int j=0; j<jMax; j++)
    {
        for (int l=0; l<lMax; l++)
        {
            write_output << uPre[l][j] << "\n";
        }
    }

    write_output.close();
    
    return 0;
}

void MatrixMultiply(double M1[20][10], double M2[10][20], double M[20][20], int order)
{
    //Multiply matrix M1 of size order x 10 and matrix M2 of size 10 x order to yield matrix M of size order x order 
    int count;
    count = 0;
    for (int i=0; i<order; i++)
    {
        for (int j=0; j<order; j++)
        {        
            for (int k=0; k<10; k++)
            {
                M[i][j] += M1[i][k]*M2[k][j];
            }
        }
    }
}

void MatrixMultiply2(double M1[20][20], double M2[20][10], double M[20][10], int order)
{
    //Multiply matrix M1 of size order x order and matrix M2 of size order x 10 to yield matrix M of size order x 10
    for (int i=0; i<order; i++)
    {
        for (int j=0; j<10; j++)
        {        
            for (int k=0; k<order; k++)
            {
                M[i][j] += M1[i][k]*M2[k][j];

            }
        }
    }
}

void MatrixInverse(double M[20][20], double M_inverse[20][20], int order)
{
    //Take inverse of Matrix M through Gaussian Jordan elimination
    double ratio;
    double M_inverse_big[20][20];
    
    // Augmenting Identity Matrix of Order n 
    for(int i=0; i<order; i++)
    {
        for(int j=0; j<order; j++)
        {
            M_inverse_big[i][j] = M[i][j];
            if (i==j)
            {
                M_inverse_big[i][j+order] = 1;
            }
            else
            {
                M_inverse_big[i][j+order] = 0;
            }
        }
    }
    // Applying Gauss Jordan Elimination 
    for(int i=0; i<order; i++)
    {
        assert(M_inverse_big[i][i] != 0.0);
        for(int j=0; j<order; j++)
        {
            if(i!=j)
            {
                ratio = M_inverse_big[j][i]/M_inverse_big[i][i];
                for(int k=0; k<2*order; k++)
                {
                    M_inverse_big[j][k] = M_inverse_big[j][k] - ratio*M_inverse_big[i][k];
                }
            }
        }
    }
    // Row Operation to Make Principal Diagonal to 1 
    for(int i=0; i<order;i++)
    {
        for(int j=order; j<2*order; j++)
        {
            M_inverse[i][j-order] = M_inverse_big[i][j]/M_inverse_big[i][i];
        }
    }

}

void LeastSquares(double uInitialize[], double xj, double dx, int order)
{
    //Compute initial condition using Least Squares method
    double x;
    double y[10]={0};
    double bigX[10][20]={0};
    double bigXT[20][10]={0};
    double bigXprod[20][20]={0};
    double bigXprodInv[20][20]={0};
    double hugeX[20][10]={0};
    for (int i=0; i<10; i++)
    {
        x = xj-dx/2.0+i*dx/9.0;
        y[i] = sin(x);
        for (int l=0; l<order; l++)
        {
            bigX[i][l] = LegendreP(l,2.0*(x-xj)/dx);
            bigXT[l][i] = LegendreP(l,2.0*(x-xj)/dx);
        }
    }
    MatrixMultiply(bigXT,bigX,bigXprod,order);
    MatrixInverse(bigXprod,bigXprodInv,order);
    MatrixMultiply2(bigXprodInv,bigXT,hugeX,order);

    for (int i=0; i<order; i++)
    {      
        for (int j=0; j<10; j++)
        {
            uInitialize[i] += hugeX[i][j]*y[j];
        }
    }
}

void MatrixMultiply3(double M1[20][20], double M2[20][20], double M[20][20], int order)
{
    //Multiply matrix M1 of size order x order and matrix M2 of size order x order to yield matrix M of size order x order
    for (int i=0; i<order; i++)
    {
        for (int j=0; j<order; j++)
        {        
            for (int k=0; k<order; k++)
            {
                M[i][j] += M1[i][k]*M2[k][j];
            }
        }
    }
}

double LegendreP(int n, double x)
{
    //Calculate P_n(x), the n-th order Legendre polynomial at x
    assert (n>=0);
    switch (n)
    {
        case 0:
            return 1;
            break;
        case 1:
            return x;
            break;
        default:
            return ((2.0*n-1.0)*x*LegendreP(n-1,x)-(n-1)*LegendreP(n-2,x))/n;
            break;
    }
}

double LegendrePorthonormal(int n, double x)
{
    //Make P_n(x), the n-th order Legendre polynomial at x, normal
    assert (n>=0);
    return sqrt((2.0*n+1.0)/2.0)*LegendreP(n,x);
}

double LegendreP_derivative(int n, double x)
{
    //Calculate P'_n(x), the first derivative of the n-th order Legendre polynomial at x
    assert (n>=0);
    switch (n)
    {
        case 0:
            return 0;
            break;
        case 1:
            return 1;
            break;
        default:
            return ((2.0*n-1.0)*(x*LegendreP_derivative(n-1,x)+LegendreP(n-1,x))-(n-1)*LegendreP_derivative(n-2,x))/n;
    }
}

double LegendrePorthonormal_derivative(int n, double x)
{
    //Make P'_n(x), the first derivative of the n-th order Legendre polynomial at x, normal
    assert(n>=0);
    return sqrt((2.0*n+1.0)/2.0)*LegendreP_derivative(n,x);
}

void roots_initial(int n, double x[])
{
    //Initial guess of roots of Legendre Polynomial of order n
    assert(n>=0);
    for (int i=0; i<n; i++)
    {
        x[i] = (1.0-1.0/(8.0*pow(n,2.0))+1.0/(8.0*pow(n,3.0)))*cos(M_PI*(4.0*(i+1.0)-1.0)/(4.0*n+2.0));
    }

}

void roots_final(int n, double x[])
{
    //Newton's method to find roots of Legendre Polynomial of order n
    double x_root;
    assert(n>=0);
    for (int i=0; i<n; i++)
    {
        x_root = x[i];
        while (fabs(LegendreP(n,x_root))>1.0e-10)
        {
            x_root = x_root-LegendreP(n,x_root)/LegendreP_derivative(n,x_root);
        }
        x[i] = x_root;
    }

}

void weights(int n, double x[], double w[])
{
    //Calculate weights for order n Gaussian Quadrature integration using Legendre points (note order of 2n+1 needed to integrate function of order n)
    assert(n>=0);
    double x_root;
    for (int i=0; i<n; i++)
    {
        x_root = x[i];
        w[i] = 2.0/((1.0-pow(x_root,2.0))*pow(LegendreP_derivative(n,x_root),2.0));
    }
}

double integrate(int n_func1, int n_func2, int n, double x[], double w[])
{
    //Integrate two Legendre polynomials of order n_func1 and n_func2
    assert(n>=0);
    double y=0;
    for (int i=0; i<n; i++)
    {
        y+=w[i]*LegendreP(n_func1,x[i])*LegendreP(n_func2,x[i]);
    }
    return y;
}

double integrateS(int n_func1, int n_func2, int n, double x[], double w[])
{
    //Integrate the derivative of a Legendre polynomial of order n_func1 with a Legendre polynomial of order n_func2
    assert(n>=0);
    double y=0;
    for (int i=0; i<n; i++)
    {
        y+=w[i]*LegendreP_derivative(n_func1,x[i])*LegendreP(n_func2,x[i]);
    }
    return y;
}


