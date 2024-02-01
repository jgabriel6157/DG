#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>

void readFile(std::string filename, std::string argName[], std::string argString[], int numberOfVariables);
int assignInt(std::string varString);
double assignDouble(std::string varString);

double getFunction(std::string basis, int n, double x, bool derivative=false);
double LegendreP(int n, double x);
double LegendreP_derivative(int n, double x);
double LegendrePorthonormal(int n, double x);
double LegendrePorthonormal_derivative(int n, double x);

void LeastSquares(double* uInitialize, double xj, double dx, std::string basis, int order);
void MatrixMultiply(double** M1, double** M2, double** M, int rowM1, int middleSize, int columnM2);
void MatrixMultiply(double** M, double* A1, double* A, int rowM1, int middleSize);
void MatrixInverse(double** M, double** M_inverse, int order);
void roots_initial(int n, double x[]);
void roots_final(int n, double x[]);
void weights(int n, double x[], double w[]);
double integrate(std::string basis, int n_func1, bool derivative1, int n_func2, bool derivative2, int n, double x[], double w[]);


int main(int argc, char* argv[])
{
    std::string argName[8] = {"a","jMax","lMax","tMax","quadratureOrder","length","dt","basis"};
    std::string argString[8];
    readFile("input.txt",argName,argString,8);

    double a = assignDouble(argString[0]);
    int jMax = assignInt(argString[1]);
    int lMax = assignInt(argString[2]);
    int tMax = assignInt(argString[3]);
    int quadratureOrder = assignInt(argString[4]);
    double length = assignDouble(argString[5]);
    double dt = assignDouble(argString[6]);
    std::string basis = argString[7];
    double dx = length/jMax;

    double** M = new double* [lMax];
    double** M_inv = new double* [lMax];
    double** S = new double* [lMax];
    double** F1 = new double* [lMax];
    double** F2 = new double* [lMax];
    double** M_invS = new double* [lMax];
    double** M_invF1 = new double* [lMax];
    double** M_invF2 = new double* [lMax];
    double** uPre = new double* [lMax];
    double** uPost = new double* [lMax];
    for (int l=0; l<lMax; l++)
    {
        M[l] = new double [lMax];
        M_inv[l] = new double [lMax];
        S[l] = new double [lMax];
        F1[l] = new double [lMax];
        F2[l] = new double [lMax];
        M_invS[l] = new double [lMax];
        M_invF1[l] = new double [lMax];
        M_invF2[l] = new double [lMax];
        uPre[l] = new double [jMax];
        uPost[l] = new double [jMax];
    }
    double* uInitialize = new double [lMax];
    double x_roots[quadratureOrder], w[quadratureOrder];
    double xj;
    
    std::ofstream write_output("Output.csv");
    assert(write_output.is_open());

    roots_initial(quadratureOrder,x_roots);
    roots_final(quadratureOrder,x_roots);
    weights(quadratureOrder,x_roots,w);

    for (int i=0; i<lMax; i++)
    {
        for (int j=0; j<lMax; j++)
        {
            M[i][j] = integrate(basis,i,false,j,false,quadratureOrder,x_roots,w)/2.0;
            S[i][j] = integrate(basis,i,true,j,false,quadratureOrder,x_roots,w);
            F1[i][j] = getFunction(basis,i,1)*getFunction(basis,j,1);
            F2[i][j] = getFunction(basis,i,-1)*getFunction(basis,j,1);
            // if (M[i][j] < 1e-10)
            // {
            //     M[i][j] = 0;
            // }
        }
    }
    
    MatrixInverse(M,M_inv,lMax);
    for (int l=0; l<lMax; l++)
    {
        delete[] M[l];
    }
    delete[] M;

    MatrixMultiply(M_inv,S,M_invS,lMax,lMax,lMax);
    MatrixMultiply(M_inv,F1,M_invF1,lMax,lMax,lMax);
    MatrixMultiply(M_inv,F2,M_invF2,lMax,lMax,lMax);

    for (int l=0; l<lMax; l++)
    {
        delete[] M_inv[l];
        delete[] S[l];
        delete[] F1[l];
        delete[] F2[l];
    }
    delete[] M_inv;
    delete[] S;
    delete[] F1;
    delete[] F2;

    for (int j=0; j<jMax; j++)
    {
        xj = j*dx+dx/2.0;
        LeastSquares(uInitialize, xj, dx, basis, lMax);
        for (int l=0; l<lMax; l++)
        {
            uPre[l][j] = uInitialize[l];
        }
    }
    delete[] uInitialize;

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
                    if (j==0)
                    {
                        uPost[l][j]+=M_invF2[l][i]*uPre[i][jMax-1];
                    }
                    else
                    {
                        uPost[l][j]+=M_invF2[l][i]*uPre[i][j-1];
                    }
                }
                uPost[l][j]*=a;
                uPost[l][j]*=dt;
                uPost[l][j]/=dx;
                uPost[l][j]+=uPre[l][j];
            }
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

    for (int l=0; l<lMax; l++)
    {
        delete[] M_invS[l];
        delete[] M_invF1[l];
        delete[] M_invF2[l];
        delete[] uPre[l];
        delete[] uPost[l];
    }
    delete[] M_invS;
    delete[] M_invF1;
    delete[] M_invF2;
    delete[] uPre;
    delete[] uPost;

    return 0;
}

void readFile(std::string filename, std::string argName[], std::string argString[], int numberOfVariables)
{
    std::ifstream inputFile(filename);
    assert(inputFile.is_open());

    std::string word,value,dum1,valueCheck;

    while (inputFile >> word)
    {
        for (int i=0;i<numberOfVariables;i++)
        {
            if (word==argName[i])
            {
                inputFile >> word >> value;
                argString[i] = value;
            }
        }
    } 

    inputFile.close();
}

int assignInt(std::string varString)
{
    int number = 1;
    double const pi = M_PI;
    std::string valueCheck = "*";
    std::string value = varString;
    std::string dum1;
    while (valueCheck.find("*") != std::string::npos)
    {
        valueCheck = value;
        dum1 = value.substr(0,value.find("*"));
        if (dum1 == "pi")
        {
            number*=pi;
        }
        else
        {
            number*=stoi(dum1);
        }
        value = value.substr(value.find("*")+1);
    }
    return number;
}

double assignDouble(std::string varString)
{
    double number = 1.0;
    double const pi = M_PI;
    std::string valueCheck = "*";
    std::string value = varString;
    std::string dum1;
    while (valueCheck.find("*") != std::string::npos)
    {
        valueCheck = value;
        dum1 = value.substr(0,value.find("*"));
        if (dum1 == "pi")
        {
            number*=pi;
        }
        else
        {
            number*=stod(dum1);
        }
        value = value.substr(value.find("*")+1);
    }
    return number;
}

void MatrixMultiply(double** M1, double** M2, double** M, int rowM1, int middleSize, int columnM2)
{
    for (int i=0; i<rowM1; i++)
    {
        for (int j=0; j<columnM2; j++)
        {
            M[i][j] = 0;        
            for (int k=0; k<middleSize; k++)
            {
                M[i][j] += M1[i][k]*M2[k][j];
            }
        }
    }
}

void MatrixMultiply(double** M, double* A1, double* A, int rowM1, int middleSize)
{
    for (int i=0; i<rowM1; i++)
    {
        A[i] = 0;
        for (int j=0; j<middleSize; j++)
        {
            A[i] += M[i][j]*A1[j];
        }
    }
}

void MatrixInverse(double** M, double** M_inverse, int order)
{
    //Take inverse of Matrix M through Gaussian Jordan elimination
    double ratio;
    double M_inverse_big[order][order*2];
    
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

void LeastSquares(double* uInitialize, double xj, double dx, std::string basis, int order)
{
    //Compute initial condition using Least Squares method
    double x;
    double* y = new double [10];
    double** bigX = new double* [10];
    double** bigXT = new double* [order];
    double** bigXprod = new double* [order];
    double** bigXprodInv = new double* [order];
    double** hugeX = new double* [order];
    for (int i=0; i<10; i++)
    {
        bigX[i] = new double [order];
    }
    for (int l=0; l<order; l++)
    {
        bigXT[l] = new double [10];
        bigXprod[l] = new double [order];
        bigXprodInv[l] = new double [order];
        hugeX[l] = new double [10];
    }

    for (int i=0; i<10; i++)
    {
        x = xj-dx/2.0+i*dx/9.0;
        y[i] = sin(x);
        for (int l=0; l<order; l++)
        {
            bigX[i][l] = getFunction(basis,l,2.0*(x-xj)/dx);
            bigXT[l][i] = getFunction(basis,l,2.0*(x-xj)/dx);
        }
    }

    MatrixMultiply(bigXT,bigX,bigXprod,order,10,order);
    for (int i=0; i<10; i++)
    {
        delete[] bigX[i];
    }
    delete[] bigX;

    MatrixInverse(bigXprod,bigXprodInv,order);
    for (int l=0; l<order; l++)
    {
        delete[] bigXprod[l];
    }
    delete[] bigXprod;

    MatrixMultiply(bigXprodInv,bigXT,hugeX,order,order,10);
    for (int l=0; l<order; l++)
    {
        delete[] bigXprodInv[l];
        delete[] bigXT[l];
    }
    delete[] bigXprodInv;
    delete[] bigXT;

    MatrixMultiply(hugeX,y,uInitialize,order,10);
    for (int l=0; l<order; l++)
    {
        delete[] hugeX[l];
    }
    delete[] hugeX;
    delete[] y;
}

double getFunction(std::string basis, int n, double x, bool derivative)
{
    if (basis=="Legendre")
    {
        if (!derivative)
        {
            return LegendreP(n,x);
        }
        else
        {
            return LegendreP_derivative(n,x);
        }
    }
    else if (basis=="LegendreOrthonormal")
    {
        if (!derivative)
        {
            return LegendrePorthonormal(n,x);
        }
        else
        {
            return LegendrePorthonormal_derivative(n,x);
        }
    }
    else
    {
        return 0;
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

double integrate(std::string basis, int n_func1, bool derivative1, int n_func2, bool derivative2, int n, double x[], double w[])
{
    //Integrate two Legendre polynomials of order n_func1 and n_func2
    assert(n>=0);
    double y=0;
    for (int i=0; i<n; i++)
    {
        y+=w[i]*getFunction(basis,n_func1,x[i],derivative1)*getFunction(basis,n_func2,x[i],derivative2);
    }
    return y;
}



