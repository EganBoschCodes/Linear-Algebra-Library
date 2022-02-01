#include "Complex.h"



Complex& Complex::operator=(const double d) {
	_real = d;
	_im = 0;
	return *this;
}
Complex& Complex::operator=(const int i) {
	_real = (double)i;
	_im = 0;
	return *this;
}
Complex& Complex::operator=(const Complex c) {
	_real = c.real();
	_im = c.im();
	return *this;
}

Complex Complex::operator+(const double d) {
	return Complex(_real + d, _im);
}
Complex Complex::operator+(const int i) {
	return Complex(_real + i, _im);
}
Complex Complex::operator+(const Complex c) {
	return Complex(_real + c.real(), _im + c.im());
}

Complex Complex::operator-() {
	return Complex(-_real, -_im);
}
Complex Complex::operator-(const double d) {
	return Complex(_real - d, _im);
}
Complex Complex::operator-(const int i) {
	return Complex(_real - i, _im);
}
Complex Complex::operator-(const Complex c) {
	return Complex(_real - c.real(), _im - c.im());
}

Complex Complex::operator*(const double d) {
	return Complex(_real * d, _im * d);
}
Complex Complex::operator*(const int i) {
	return Complex(_real * i, _im * i);
}
Complex Complex::operator*(const Complex c) {
	return Complex(_real * c.real() - _im * c.im(), _real * c.im() + _im * c.real());
}

Complex Complex::operator/(const double d) {
	return Complex(_real / d, _im / d);
}
Complex Complex::operator/(const int i) {
	return Complex(_real / i, _im / i);
}
Complex Complex::operator/(const Complex c) {
	Complex cc = *this * c.conj();
	cc = cc / (c.real() * c.real() - c.im() * c.im());
	return cc;
}



Complex& Complex::operator+=(const double d) {
	_real += d;
	return *this;
}
Complex& Complex::operator+=(const int d) {
	_real += d;
	return *this;
}
Complex& Complex::operator+=(const Complex c) {
	_real += c.real();
	_im += c.im();
	return *this;
}
Complex& Complex::operator-=(const double d) {
	_real -= d;
	return *this;
}
Complex& Complex::operator-=(const int d) {
	_real -= d;
	return *this;
}
Complex& Complex::operator-=(const Complex c) {
	_real -= c.real();
	_im -= c.im();
	return *this;
}
Complex& Complex::operator*=(const double d) {
	_real *= d;
	_im *= d;
	return *this;
}
Complex& Complex::operator*=(const int d) {
	_real *= d;
	_im *= d;
	return *this;
}
Complex& Complex::operator*=(const Complex c) {
	*this = *this * c;
	return *this;
}
Complex& Complex::operator/=(const double d) {
	_real /= d;
	_im /= d;
	return *this;
}
Complex& Complex::operator/=(const int d) {
	_real /= d;
	_im /= d;
	return *this;
}
Complex& Complex::operator/=(const Complex c) {
	*this = *this / c;
	return *this;
}


bool Complex::operator==(const Complex c) {
	return soft_round(_real) == soft_round(c.real()) && soft_round(_im) == soft_round(c.im());
}
bool Complex::operator!=(const Complex c) {
	return !(*this == c);
}
bool Complex::operator<(const Complex c) {
	return _real < c.real();
}
bool Complex::operator>(const Complex c) {
	return _real > c.real();
}
bool Complex::operator<=(const Complex c) {
	return _real <= c.real();
}
bool Complex::operator>=(const Complex c) {
	return _real <= c.real();
}



Complex Complex::pow(double d) {
	double newMag = std::pow(mag(), d);
	double newAngle = angle() * d;

	double newR = cos(newAngle) * newMag;
	double newI = sin(newAngle) * newMag;
	//if (newR < 0) {
		//newR *= -1;
		//newI *= -1;
	//}
	return Complex(newR, newI);
}
Complex Complex::sqrt() {
	return pow(0.5);
}


std::string stringify(double value) {
	std::ostringstream os;
	os << value;
	return os.str();
}

double soft_round(double value) {
	return round(value * 10000000) / 10000000;
}