#ifndef VECTORHEADERDEF
#define VECTORHEADERDEF

class Vector
{
private:
double* mData; // data stored in vector
int mSize; // size of vector
public:
Vector(const Vector& otherVector);
Vector(int size);
~Vector();
int GetSize() const;
double& operator[](int i); // zero-based indexing
// read-only zero-based indexing
double Read(int i) const;
Vector& operator=(const Vector& otherVector);
Vector operator+(double constant) const; // unary +
Vector operator-(double constant) const; // unary -
Vector operator+(const Vector& v1) const; // binary +
Vector operator-(const Vector& v1) const; // binary -
// scalar multiplication
Vector operator*(double a) const;
// vector multiplication (dot product)
double operator*(const Vector& v1) const;
// p-norm method
double CalculateNorm(int p=2) const;
// declare length function as a friend
friend int length(const Vector& v);
};

// Prototype signature of length() friend function
int length(const Vector& v);

#endif
