


#ifndef COMPLEX
#define COMPLEX

#include <math.h>
#include <iostream>
#include <string.h>
#include <sstream>



using namespace std; 

std::string stringify(double value);
double soft_round(double value);

class Complex {
public:
	Complex() : _real(0), _im(0) {}
	Complex(double d) : _real(d), _im(0) {}
	Complex(double d, double i) : _real(d), _im(i) {}

	double real() const { return _real; }
	void setReal(double r) { _real = r; }
	double im() const { return _im; }
	void setIm(double i) { _im = i; }
	double mag() const { return std::sqrt(_real * _real + _im * _im); }
	double angle() const { return atan2(_im, _real); }
	Complex conj() const { return Complex(_real, -_im); }


	Complex& operator=(const double d);
	Complex& operator=(const int i);
	Complex& operator=(const Complex c);
	 
	Complex operator+(const double d);
	Complex operator+(const int d);
	Complex operator+(const Complex d);
	Complex operator-();
	Complex operator-(const double d);
	Complex operator-(const int d);
	Complex operator-(const Complex d);
	Complex operator*(const double d);
	Complex operator*(const int d);
	Complex operator*(const Complex d);
	Complex operator/(const double d);
	Complex operator/(const int d);
	Complex operator/(const Complex d);

	Complex& operator+=(const double d);
	Complex& operator+=(const int d);
	Complex& operator+=(const Complex d);
	Complex& operator-=(const double d);
	Complex& operator-=(const int d);
	Complex& operator-=(const Complex d);
	Complex& operator*=(const double d);
	Complex& operator*=(const int d);
	Complex& operator*=(const Complex d);
	Complex& operator/=(const double d);
	Complex& operator/=(const int d);
	Complex& operator/=(const Complex d);

	bool operator==(const Complex d);
	bool operator!=(const Complex d);
	bool operator<(const Complex d);
	bool operator>(const Complex d);
	bool operator<=(const Complex d);
	bool operator>=(const Complex d);

	Complex sqrt();
	Complex pow(double pow);
	Complex round() { return Complex(std::round(real()), std::round(im())); }
	Complex softRound() { return Complex(soft_round(real()), soft_round(im())); }

	void print() { cout << "(" << _real << ", " << _im << ")" << endl; }
	string toString() { return _im == 0 ? (stringify(soft_round(_real))) : "(" + stringify(soft_round(_real)) + ", " + stringify(soft_round(_im)) + ")"; }

private:
	double _real;
	double _im;
};

const Complex I(0, 1);

#endif