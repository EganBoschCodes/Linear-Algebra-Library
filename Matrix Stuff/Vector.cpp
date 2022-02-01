#include "Vector.h"
#include "Matrix.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <random>
#include <time.h>
#include "Complex.h"



using namespace std;

Complex fRand(Complex fMin, Complex fMax)
{
	Complex f = (Complex)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}


Vector::Vector(std::string s) {
	vector<string> spl = split(s, ' ');
	for (int i = 0; i < spl.size(); i++) {
		data.push_back(stod(spl[i]));
	}
}

Vector::Vector(Complex d, Complex d2) {
	data.push_back(d);
	data.push_back(d2);
}

Vector::Vector(Complex d, Complex d2, Complex d3) {
	data.push_back(d);
	data.push_back(d2);
	data.push_back(d3);
}



Complex Vector::dot(Vector v) {
	if (v.size() != size()) {
		return 0;
	}
	Complex sum = 0;
	for (int i = 0; i < size(); i++) {
		sum += v[i] * data[i];
	}
	return sum;
}

Vector Vector::cross(Vector v) {
	if (v.size() != size() || v.size() != 3) {
		return Vector();
	}
	Vector n = *this;
	Vector k;
	k.push(n[1] * v[2] - n[2] * v[1]);
	k.push(n[2] * v[0] - n[0] * v[2]);
	k.push(n[0] * v[1] - n[1] * v[0]);
	return k;
}

void Vector::trimTo(int i) {
	while (size() > i && i >= 0) {
		data.pop_back();
	}
}

const string Vector::toString() {
	stringstream s;
	for (int i = 0; i < size(); i++) { s << data[i].real() << " "; }
	s << endl << "Magnitude: " << (!(*this)).real();
	if (size() == 2) {
		double deg = atan2(data[1].real(), data[0].real()) * 180 / 3.141592653;
		s << ", Degree: " << (deg < 0 ? deg+360 : deg);
	}
	s << endl;
	return s.str();
}

void const Vector::print() { cout << toString() << endl; }

Complex& Vector::operator[](int i) { return data[i%data.size()]; }

Vector& Vector::operator=(const Vector& v) {
	data = v.data;
	return *this;
}

Vector Vector::operator+(const Vector& shift) {
	if (shift.size() != size()) {
		return *this;
	}
	Vector v(*this);
	for (int i = 0; i < size(); i++) {
		v.data[i] += Vector(shift).get(i);
	}
	return v;
}

void Vector::operator+=(const Vector& shift) {
	for (int i = 0; i < size(); i++) {
		data[i] += shift.get(i);
	}
}

void Vector::operator-=(const Vector& shift) {
	for (int i = 0; i < size(); i++) {
		data[i] -= shift.get(i);
	}
}

Vector& Vector::operator+(Complex newVal) {
	data.push_back(newVal);
	return *this;
}

Vector Vector::operator*(Complex scale) {

	//cout << scale << endl << endl;
	Vector v(*this);
	for (int i = 0; i < size(); i++) {
		v.data[i] *= scale;
	}
	return v;
}

Vector Vector::operator*(double scale) {

	//cout << scale << endl << endl;
	Vector v(*this);
	for (int i = 0; i < size(); i++) {
		v.data[i] *= scale;
	}
	return v;
}

Vector Vector::operator/(Complex scale) {
	Vector v(*this);
	for (int i = 0; i < size(); i++) {
		v.data[i] /= scale;
	}
	return v;
}

Complex Vector::operator*(const Vector& scale) {
	return dot(scale);
}

Vector Vector::operator-(const Vector& shift) {
	if (shift.size() != size()) {
		return *this;
	}
	Vector v(*this);
	for (int i = 0; i < size(); i++) {
		v.data[i] -= shift.get(i);
	}
	return v;
}

Vector Vector::operator-() {
	Vector v = *this;
	for (int i = 0; i < size(); i++) {
		v.data[i] *= -1;
	}
	return v;
}

Vector Vector::operator~() {
	Vector v = *this;
	Complex mag = !(*this);
	for (int i = 0; i < size(); i++) {
		v.data[i] /= mag;
	}
	return v;
}

Complex Vector::operator!() {
	double sum = 0;
	for (int i = 0; i < size(); i++) {
		sum += data[i].mag() * data[i].mag();
	}
	return sqrt(sum);
}

Matrix Vector::operator&&(const Vector& v) {
	Matrix m(0, v.size());
	return (m && *this && v);

}

Vector Vector::operator||(Vector vec) {
	Vector thisNorm = ~(*this);
	return thisNorm * (vec * thisNorm);
}

bool Vector::operator==(const Vector& comp) {
	if (comp.size() != size()) { return false; }
	print();
	for (int i = 0; i < size(); i++) {
		cout << comp.get(i).real() << " ";
	}
	for (int i = 0; i < size(); i++) {
		if (abs(comp.get(i).real() - get(i).real()) < 0.000001) {
			return false;
		}
	}
	return true;
}


int numWholes(Vector& v) {
	int a = 0;
	for (int i = 0; i < v.size(); i++) {
		if (abs(round(v[i].real()) - v[i].real()) < 0.000001) {
			v[i].setReal(round(v[i].real()));
			a += 1;
		}
	}
	return a;
}

void Vector::wholeNumberize() {
	int multiplier = 2;
	int wholes = numWholes(*this);
	while (wholes < size() && multiplier < 1000) {
		*this = *this * multiplier;
		int newwholes = numWholes(*this);
		if (newwholes > wholes) {
			//cout << "New Multiple: " << multiplier << endl;
			wholes = newwholes;
			//cout << (size() - wholes) << " more decimals." << endl;
			//print();
			multiplier = 2;
		}
		else {
			*this = *this / multiplier;
			multiplier++;
		}
	}

	return;
}

bool Vector::isZeroVector() {
	bool a = true;
	for (int i = 0; i < size(); i++) {
		a = a && (get(i).mag() < 0.0000001);
	}
	return a;
}

int Vector::firstNonZeroEntry() {
	int i = -1;
	Complex val = 0;
	while (i + 1 < size() && val == 0) {
		val = data[i + 1];
		i++;
	}
	return i + (val == 0 ? 1 : 0);
}

void Vector::deleteEntry(int n) {
	data.erase(data.begin() + n);
}







Vector randomVector(int dim) {
	Vector v(dim);
	for (int i = 0; i < dim; i++) {
		v[i] = fRand(-1, 1);
	}
	return v;
}

Vector VectorFromRadAndMag(double angle, double mag) {
	return Vector(cos(angle) * mag,sin(angle) * mag);
}

Vector VectorFromDegreeAndMag(double angle, double mag) {
	return VectorFromRadAndMag(angle / 180 * 3.141592653, mag);
}