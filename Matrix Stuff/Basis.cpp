#include "Basis.h"
#include "Vector.h"
#include "Matrix.h"

using namespace std;

Basis::Basis() : dim(0), verified(false), orthonormal(false){}

Basis::Basis(Matrix m) : dim(m.numCols()), verified(false), orthonormal(false) {
	for (int i = 0; i < dim; i++) {
		data.push_back(Vector(m.getColumn(i)));
	}

	removeRepeats();
}

Basis::Basis(const Basis& b) : dim(b.dim), data(b.data), verified(false), orthonormal(false) {}

Basis::Basis(string s): dim(1), verified(false), orthonormal(false) {
	vector<string> spl = split(s, ','); 
	data.push_back(Vector(spl[0]));
	for (int i = 1; i < spl.size(); i++) {
		*this += Vector(spl[i]);
	}
}



Basis& Basis::operator=(const Basis& b) {
	data = b.data;
	dim = b.dim;
	return *this;
}



Basis& Basis::operator+(Vector v) {
	if (!(*this && v) > 0) {
		dim++;
		data.push_back(v);
		verified = false;
		orthonormal = false;
	}
	return *this;
}

void Basis::operator+=(Vector v) {
	*this = *this + v;
}

Basis& Basis::operator+(Basis b) {
	for (int i = 0; i < b.dim; i++) {
		*this += b.getVector(i);
	}

	removeRepeats();
	return *this;
}

void Basis::operator+=(Basis b) {
	*this = *this + b;
}



Vector Basis::operator&&(Vector v) {
	Vector per = *this || v;
	return  (v - per);
}

Vector Basis::operator||(Vector v) {
	if (dim == 0) {
		return Vector(v.size());
	}


	Basis b = isOrthoNormal() ? getOrthoNormal() : *this;

	Matrix m(b);
	Matrix mt(~m);
	return m * mt * v;
}


Basis Basis::operator!() {

	if (dim == 0) {
		return Basis();
	}

	if (data[0].size() == 0) {
		return Basis();
	}

	int size = data[0].size();
	if (size == dim) {
		return Basis();
	}

	srand(data[0][0].real());
	Basis copy(*this);

	Vector v = randomVector(size);
	while (!(copy && v) == 0) {
		v = randomVector(size);
	}
	copy += ~(copy && v);

	for (int i = 0; i < size - dim - 1; i++) {
		Vector v = randomVector(size);
		while (!(copy && v) == 0) {
			v = randomVector(size);
		}
		copy += ~(copy && v);
	}

	copy.verified = false;
	copy.orthoNormalize();

	Basis b;

	for (int i = 0; i < size - dim; i++) {
		b += copy.getVector(i + dim);
	}

	return b;
}

Basis Basis::getOrthoNormal() {
	if (dim == 0) {
		return *this;
	}

	if (isOrthoNormal()) {
		return *this;
	}

	Basis temp;
	temp += ~(data[0]);
	for (int i = 1; i < dim; i++) {
		Vector v = temp && data[i];
		if (!v.isZeroVector()) {
			temp += ~v;
		}
	}
	return temp;
}

void Basis::orthoNormalize() {
	*this = getOrthoNormal();
}


bool Basis::isOrthoNormal() {
	if (verified) {
		return orthonormal;
	}
	else {
		return isOrthoNormal_private();
	}
}

bool Basis::isOrthoNormal_private() {
	if (data.size() == 0) {
		return true;
	}

	orthonormal = true;
	for (int i = 0; i < data.size(); i++) {
		orthonormal = orthonormal && (abs(1 - (!data[i]).mag()) < 0.00001);
		for (int b = 0; b < i; b++) {
			orthonormal = orthonormal && ((data[i] * data[b]).mag() < 0.00001);
		}
	}

	verified = true;

	return orthonormal;
}



void Basis::print() {
	cout << "Dimensions: " << dim << endl;
	cout << "-----------------" << endl;
	for (int i = 0; i < data.size(); i++) {
		data[i].print();
	}
	cout << endl;
}

Vector Basis::getVector(int i) { return data[i]; }

bool Basis::fillsDimension() { return data.size() > 0 ? dim == data[0].size() : false; }




void Basis::removeRepeats() {
	if (dim == 0 || (verified && orthonormal)) {
		return;
	}


	Basis temp;
	temp += data[0];
	for (int i = 1; i < dim; i++) {
		Vector v = temp && data[i];
		if (!v.isZeroVector()) {
			temp += data[i];
		}
	}

	*this = temp;
}