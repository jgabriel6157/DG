#include <iostream>
#include <cmath>
#include <cassert>
#include "Vector.hxx"

// Overridden copy constructor
// Allocates memory for new vector, and copies
// entries of other vector into it
Vector::Vector(const Vector& otherVector)
{
    mSize = otherVector.GetSize();
    mData = new double [mSize];
    for (int i=0; i<mSize; i++)
    {
        mData[i] = otherVector.mData[i];
    }
}

// Constructor for vector of a given size
// Allocates memory, and initialises entries
// to zero
Vector::Vector(int size)
{
    assert(size > 0);
    mSize = size;
    mData = new double [mSize];
    for (int i=0; i<mSize; i++)
    {
        mData[i] = 0.0;
    }
}

// Overridden destructor to correctly free memory
Vector::~Vector() 
{
    delete[] mData;
}

// Method to get the size of a vector
int Vector::GetSize() const
{
    return mSize;
}

void Vector::Print()
{
    for (int i=0; i<mSize; i++)
    {
        std::cout << mData[i] << "\n";
    }
}

// Overloading square brackets
// Note that this uses ‘zero-based’ indexing,
// and a check on the validity of the index
double& Vector::operator[](int i)
{
    assert(i > -1);
    assert(i < mSize);
    return mData[i];
}

// Read-only variant of []
// Note that this uses ‘zero-based’ indexing,
// and a check on the validity of the index
double Vector::Read(int i) const
{
    assert(i > -1);
    assert(i < mSize);
    return mData[i];
}

// Overloading the assignment operator
Vector& Vector::operator=(const Vector& otherVector)
{
    assert(mSize == otherVector.mSize);
    for (int i=0; i<mSize; i++)
    {
        mData[i] = otherVector.mData[i];
    }
    return *this;
}

// Overloading the unary + operator
Vector Vector::operator+(double constant) const
{
    Vector v(mSize);
    for (int i=0; i<mSize; i++)
    {
        v[i] = constant+mData[i];
    }
    return v;
}

// Overloading the unary - operator
Vector Vector::operator-(double constant) const
{
    Vector v(mSize);
    for (int i=0; i<mSize; i++)
    {
        v[i] = mData[i]-constant;
    }
    return v;
}

// Overloading the binary + operator
Vector Vector::operator+(const Vector& v1) const
{
    assert(mSize == v1.mSize);
    Vector v(mSize);
    for (int i=0; i<mSize; i++)
    {
        v[i] = mData[i] + v1.mData[i];
    }
    return v;
}

// Overloading the binary - operator
Vector Vector::operator-(const Vector& v1) const
{
    assert(mSize == v1.mSize);
    Vector v(mSize);
    for (int i=0; i<mSize; i++)
    {
        v[i] = mData[i] - v1.mData[i];
    }
    return v;
}

// Overloading scalar multiplication
Vector Vector::operator*(double a) const
{
    Vector v(mSize);
    for (int i=0; i<mSize; i++)
    {
        v[i] = a*mData[i];
    }
    return v;
}

// Overloading vector multiplication
double Vector::operator*(const Vector& v1) const
{
    assert(mSize == v1.mSize);
    double result=0;
    for (int i=0; i<mSize; i++)
    {
        result += mData[i] * v1.mData[i];
    }
    return result;
}

// Method to calculate norm (with default value p=2)
// corresponding to the Euclidean norm
double Vector::CalculateNorm(int p) const
{
    double norm_val, sum = 0.0;
    for (int i=0; i<mSize; i++)
    {
        sum += pow(fabs(mData[i]), p);
    }
    norm_val = pow(sum, 1.0/((double)(p)));
    return norm_val;
}

// MATLAB style friend to get the size of a vector
int length(const Vector& v)
{
    return v.mSize;
}
