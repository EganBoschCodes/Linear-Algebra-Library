#ifndef MATRIX
#define MATRIX

#include <vector>
#include <string>
#include "Complex.h"

using namespace std;

class Vector;
class Basis;

class EigenVector {
public:
	EigenVector();
	EigenVector(Complex value, Vector vector);
	~EigenVector() {}

	void print();
	Vector vector();
	Complex value();
private:
	Complex val;
	Vector* vec;
};

class Matrix {
public:
	Matrix();

	//Set Size and Fill with Zeroes
	Matrix(unsigned int r, unsigned int c);

	//Set Size and Fill with Function taking coordinate as arguments.
	Matrix(unsigned int r, unsigned int c, Complex (*f)(int, int));

	//Duplicators
	Matrix(const Matrix& old);
	Matrix(const Matrix* old);

	//Simple String Construction
	Matrix(string construc);
	~Matrix() {}

	//Construct a matrix with columns equal to the vectors in b
	Matrix(Basis b);

	//Useful Functions
	bool setRow(Complex* row, unsigned int row_num);
	bool setRow(vector<Complex> row, unsigned int row_num);
	bool setRow(Vector row, unsigned int row_num);
	bool multRow(Complex scale, unsigned int row_num);
	bool addMultRowTo(Complex scale, unsigned int src_num, unsigned int target_num);

	Vector getRow(int i) const;
	Vector getColumn(int i) const; 
	void setRow(Vector v, int i);
	void setColumn(Vector v, int i);
	void delRow(int i);
	void delColumn(int i);

	//Reduced Row Echelon Form
	bool toRREF();


	//Standard Math
	Matrix operator+(const Matrix& b);
	Matrix operator-(const Matrix& b);
	Matrix operator-();
	Matrix operator*(const Matrix& b);
	Vector operator*(Vector b);
	Matrix operator*(Complex b);
	bool operator==(const Matrix& m);

	//This appends a vector (or many vectors) as a column
	Matrix& operator||(Vector col);
	Matrix& operator||(Matrix m);

	//This appends a vector (or many vectors) as a row
	Matrix& operator&&(Vector row);
	Matrix& operator&&(Matrix m);
	
	//Inverse Matrix
	Matrix operator!();

	//Transpose
	Matrix operator~();


	//Various Matrix math
	unsigned int rank();
	bool invertible();
	Vector betaCoords(Vector v);
	Matrix betaMatrix(Matrix t);
	Matrix betaMatrix(Basis b);
	bool isOrthogonal();
	bool isSymmetric();
	bool isLinearlyIndependent();
	bool getQR(Matrix& q, Matrix& r);
	Vector getLSS(Vector b);
	Matrix getOrthoProj();
	Complex determinant();
	bool isPerpendicular(Vector b);
	Basis getImage();
	Basis getNullImage();
	Complex getVolume();
	bool swapRows(unsigned int a, unsigned int b);
	Matrix getCofactor(unsigned int r, unsigned int c);
	Matrix getAdjoint();
	bool isIdentity();
	Complex trace();
	Matrix exp();
	Matrix pow(int i);
	//Returns All EigenValues between Min and Max
	vector<Complex> getEigenValues();
	vector<EigenVector> getEigenVectors();

	//Element GET/SET
	void set(Complex d, unsigned int row, unsigned int col);
	Complex get(int r, int c) const;

	//Sizing
	int numRows() const { return n_rows; }
	int numCols() const { return n_cols; }
	int numElements() const { return n_rows * n_cols; }

	//Useful for debugging
	void print();
	void activateSteps() { shouldPrintSteps = true; }

private:
	unsigned int getLowestNonZeroRow();
	unsigned int n_rows, n_cols;
	vector<Vector> data;
	bool shouldPrintSteps = false;
	Complex rightDiag(int start);
	Complex leftDiag(int start);
};

Matrix identityMatrix(unsigned int n);

Complex CramersRule(Matrix m, Vector sol, unsigned int var);
Vector CramersRule(Matrix m, Vector sol);

Vector findEigenvector(Matrix m);


#endif