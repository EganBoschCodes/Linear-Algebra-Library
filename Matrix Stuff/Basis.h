#ifndef BASIS
#define BASIS

#include <vector>
#include <string>

using namespace std;


class Vector;
class Matrix;
class Complex;

class Basis {
public:
	Basis();
	Basis(Matrix m);
	Basis(const Basis& b);
	Basis(string s);

	Basis& operator=(const Basis& b);

	//Add Vector to Basis
	Basis& operator+(Vector v);
	void operator+=(Vector v);

	//Combine two different basis
	Basis& operator+(Basis b);
	void operator+=(Basis b);

	//Perpendicular
	Vector operator&&(Vector v);

	//Projection
	Vector operator||(Vector v);

	//Null Space
	Basis operator!();

	Basis getOrthoNormal();
	void orthoNormalize();
	bool isOrthoNormal();
	
	void print();

	int dimensions() { return dim; }
	bool fillsDimension();
	Vector getVector(int i);

private:
	int dim;
	void removeRepeats(); 
	bool isOrthoNormal_private();
	bool verified, orthonormal;
	vector<Vector> data;
};

#endif