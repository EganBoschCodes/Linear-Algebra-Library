#ifndef VECTOR
#define VECTOR

#include <vector>
#include <string>
#include <cstdarg>
#include "EquationParser.h"
#include "Complex.h"

class Matrix;

class Vector {
public:
	Vector() {}
	Vector(const Vector& old) { data = old.data; }
	Vector(std::vector<Complex> d) { data = d; }
	Vector(Vector* v) { data = v->data; delete v; }
	Vector(unsigned int d) { for (int i = 0; i < d; i++) { push(Complex()); } }
	Vector(std::string s);
	Vector(Complex d, Complex d2);
	Vector(Complex d, Complex d2, Complex d3);

	void push(Complex d) { data.push_back(d); }
	void insert(Complex d, int index) { data.insert(data.begin() + index, d); }
	void trimTo(int i);
	int size() const { return data.size(); }
	std::vector<Complex> getData() const { return data; }
	Complex get(int i) const { return data[i]; }
	const std::string toString();
	const void print();


	Complex& operator[] (const int i);
	Complex get(int i) { return data[i]; }

	Vector& operator=(const Vector& old);
	Vector operator+(const Vector& shift);
	Vector& operator+(Complex newVal);
	void operator+=(const Vector& shift);
	void operator-=(const Vector& shift);
	Vector operator-(const Vector& shift);
	Vector operator*(Complex scale);
	Vector operator*(double scale);
	Vector operator/(Complex scale);
	Complex operator*(const Vector& dotter);
	Vector operator-();
	bool operator==(const Vector& comp);

	//Normalize
	Vector operator~();

	//Get Magnitude
	Complex operator!();
	Matrix operator&&(const Vector& v);

	//Get Projection onto
	Vector operator||(Vector v);

	Complex dot(Vector v);
	Vector cross(Vector v);
	void wholeNumberize();
	bool isZeroVector();
	void deleteEntry(int i);
	int firstNonZeroEntry();
	
private:
	std::vector<Complex> data;
};

bool isZeroVector(Vector v); 

Vector randomVector(int dim);

Vector VectorFromRadAndMag(Complex angle, Complex mag);

Vector VectorFromDegreeAndMag(Complex angle, Complex mag);


#endif