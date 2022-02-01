#include "Matrix.h"
#include "Vector.h"
#include <cmath>
#include "Basis.h"
#include "Complex.h"

using namespace std;



EigenVector::EigenVector() {
	val = 0;
	vec = new Vector();
}

EigenVector::EigenVector(Complex value, Vector vector) {
	val = value;
	vec = new Vector(vector);
}

void EigenVector::print() {
	cout << "Eigen-Value: " << val.toString() << endl;
	cout << "Eigen-Vector: [";
	for (int i = 0; i < vec->size(); i++) {
		cout << vec->get(i).toString();
		if (i < vec->size() - 1) { cout << " "; }
	}
	cout << "]" << endl << endl;
}

Vector EigenVector::vector() { return Vector(*vec); }

Complex EigenVector::value() { return val; }


Matrix::Matrix() : n_rows(0), n_cols(0) {};

Matrix::Matrix(unsigned int r, unsigned int c) : n_rows(r), n_cols(c) {
	for (int i = 0; i < r; i++) {
		data.push_back(vector<Complex>());
		for (int b = 0; b < c; b++) {
			data[i].push(0);
		}
	}
}

Matrix::Matrix(unsigned int r, unsigned int c, Complex (*f)(int, int)) : n_rows(r), n_cols(c) {
	for (int i = 0; i < r; i++) {
		data.push_back(vector<Complex>());
		for (int b = 0; b < c; b++) {
			data[i].push((*f)(i, b));
		}
	}
}



Matrix::Matrix(const Matrix& old): n_rows(old.n_rows), n_cols(old.n_cols) {
	data = old.data;
}

Matrix::Matrix(const Matrix* old) : n_rows(old->n_rows), n_cols(old->n_cols) {
	data = old->data;
	delete old;
}



Matrix::Matrix(string s) {
	vector<string> spl = split(s, ',');
	for (int i = 0; i < spl.size(); i++) {
		data.push_back(Vector(spl[i]));
	}
	n_cols = s.length() > 0 ? data[0].size() : 0;
	n_rows = s.length() > 0 ? spl.size() : 0;
}



Matrix::Matrix(Basis b) : n_rows(b.dimensions() != 0 ? b.getVector(0).size() : 0), n_cols(0) {
	for (int i = 0; i < n_rows; i++) {
		data.push_back(Vector());
	}
	for (int i = 0; i < b.dimensions(); i++) {
		*this = *this || b.getVector(i);
	}
}





bool Matrix::setRow(Complex* row, unsigned int row_num) {
	if (sizeof(row) != n_cols || row_num >= n_rows) {
		return false;
	}
	for (int i = 0; i < n_cols; i++) {
		data[row_num][i] = row[i];
	}
	return true;
}

bool Matrix::setRow(vector<Complex> row, unsigned int row_num) {
	if (row.size() != n_cols || row_num >= n_rows) {
		return false;
	}
	for (int i = 0; i < n_cols; i++) {
		data[row_num][i] = row[i];
	}
	return true;
}

bool Matrix::setRow(Vector row, unsigned int row_num) {
	if (row.size() != n_cols || row_num >= n_rows) {
		return false;
	}
	data[row_num] = row;
	return true;
}


void Matrix::delRow(int n) {
	data.erase(std::next(data.begin(), n), std::next(data.begin(), n + 1));
	n_rows--;
}

void Matrix::delColumn(int n) {
	for (int i = 0; i < n_rows; i++) {
		data[i].deleteEntry(n);
	}
	n_cols--;
}



bool Matrix::multRow(Complex scale, unsigned int row_num) {
	if (row_num >= n_rows) {
		return false;
	}
	for (int i = 0; i < n_cols; i++) {
		data[row_num][i] *= scale;
	}
	return true;
}

bool Matrix::addMultRowTo(Complex scale, unsigned int src_num, unsigned int target_num) {
	if (src_num >= n_rows || target_num >= n_rows) {
		return false;
	}
	for (int i = 0; i < n_cols; i++) {
		data[target_num][i] += scale * data[src_num][i];
		if (data[target_num][i].mag() < 0.0000001) { data[target_num][i] = 0; }
	}
	return true;
}



void Matrix::set(Complex d, unsigned int row, unsigned int col) { data[row][col] = d; }

Complex Matrix::get(int r, int c) const { return data[r].get(c); }





bool Matrix::toRREF() {
	int offset = 0;
	for (int i = 0; i < n_rows; i++) {
		while ((i + offset) < n_cols && data[i][i + offset] == 0) { offset++; }
		if (i + offset == n_cols && getLowestNonZeroRow() > i) {
			swapRows(i, getLowestNonZeroRow());
			offset = 0;
			i = -1;
			continue;
		}
		if (data[i][i + offset] != 0) {
			multRow(Complex(1) / data[i][i + offset], i);
			for (int b = 0; b < n_rows; b++) {
				if (b != i) {
					addMultRowTo(-data[b][i + offset] / data[i][i + offset], i, b);
				}
			}
		}
		else {
			break;
		}
	}

	for (int i = 0; i < n_rows; i++) {
		for (int b = 0; b < n_cols; b++) {
			if (data[i][b] == 0) { data[i][b] = 0; }
		}
	}

	for (int i = 1; i < n_rows; i++) {
		//cout << i << " " << getRow(i).firstNonZeroEntry() << endl;
		if (getRow(i).firstNonZeroEntry() < getRow(i - 1).firstNonZeroEntry()) {
			swapRows(i, i - 1);
			return toRREF();
		}
	}


	return true;
}




unsigned int Matrix::rank() {
	Matrix m = *this;
	m.toRREF();
	unsigned int r = 0, c = 0;
	while (r < n_rows && c < n_cols) {
		if (m.data[r][c] != 0) {
			r += 1;
		}
		c += 1;
	}
	return r;
}

bool Matrix::invertible() {
	return n_rows > 0 && n_rows == n_cols && rank() == n_rows;
}

Vector Matrix::betaCoords(Vector v) {
	Matrix m = *this;
	m || v;
	m.toRREF();
	return m.getColumn(m.n_cols - 1);
}

Matrix Matrix::betaMatrix(Matrix T) {
	Matrix B(T.n_rows, 0);
	for (int i = 0; i < n_cols; i++) {
		B = (B || betaCoords(T * getColumn(i)));
	}
	return B;
}

Matrix Matrix::betaMatrix(Basis b) {
	return betaMatrix(Matrix(b));
}




Matrix Matrix::operator+(const Matrix& b) {
	Matrix m(*this);
	if (n_cols == b.n_cols && n_rows == b.n_rows) {
		for (int r = 0; r < n_rows; r++) {
			for (int c = 0; c < n_cols; c++) {
				m.data[r][c] += b.data[r].get(c);
			}
		}
	}
	return m;
}

Matrix Matrix::operator-(const Matrix& b) {
	Matrix m(*this);
	if (n_cols == b.n_cols && n_rows == b.n_rows) {
		for (int r = 0; r < n_rows; r++) {
			for (int c = 0; c < n_cols; c++) {
				m.data[r][c] -= b.data[r].get(c);
			}
		}
	}
	return m;
}

Matrix Matrix::operator-() {
	Matrix m(*this);
	for (int r = 0; r < n_rows; r++) {
		for (int c = 0; c < n_cols; c++) {
			m.data[r][c] *= -1;
		}
	}
	return m;
}

Matrix Matrix::operator*(const Matrix& b) {
	Matrix m(n_rows, b.n_cols);
	if (b.n_rows == n_cols) {
		for (int r = 0; r < n_rows; r++) {
			for (int c = 0; c < b.n_cols; c++) {
				Complex val = getRow(r) * b.getColumn(c);
				if (val.real() - round(val.real()) < 0.0000001) { val.setReal(round(val.real())); }
				if (val.im() - round(val.im()) < 0.0000001) { val.setIm(round(val.im())); }
				m.set(val, r, c);
			}
		}
	}
	return m;
}

Vector Matrix::operator*(Vector b) {
	Vector v;
	for (int i = 0; i < n_rows; i++) {
		v.push(b * getRow(i));
	}
	return v;
}

Matrix Matrix::operator*(Complex c) {
	for (int i = 0; i < n_rows; i++) {
		for (int b = 0; b < n_cols; b++) {
			data[i][b] *= c;
		}
	}
	return *this;
}

bool Matrix::operator==(const Matrix& m) {
	if (n_rows != m.n_rows || n_cols != m.n_cols) {
		return false;
	}
	bool equivalent = true;
	for (int r = 0; r < n_rows; r++) {
		for (int c = 0; c < n_cols; c++) {
			equivalent = (get(r,c) == m.get(r,c));
			if (!equivalent) { break; }
		}
		if (!equivalent) { break; }
	}
	return equivalent;
}



Matrix& Matrix::operator&&(Vector row) {
	if (row.size() == n_cols) {
		data.push_back(row);
		n_rows = data.size();
	}
	return *this;
}

Matrix& Matrix::operator&&(Matrix m) {
	for (int i = 0; i < m.n_cols; i++) {
		*this = *this && m.getRow(i);
	}
	return *this;
}

Matrix& Matrix::operator||(Vector col) {
	if (col.size() == n_rows) {
		for (int i = 0; i < n_rows; i++) {
			data[i].push(col[i]);
			n_cols = data[i].size();
		}
	}
	
	return *this;
}

Matrix& Matrix::operator||(Matrix m) {
	for (int i = 0; i < m.n_cols; i++) {
		*this = *this || m.getColumn(i);
	}
	return *this;
}



Matrix Matrix::operator!() {
	Matrix m(*this);
	for (int i = 0; i < n_rows; i++) {
		Vector v(n_rows);
		v[i] = 1;
		m = m || v;
	}

	m.toRREF();

	Matrix inv(n_rows, 0);
	for (int i = 0; i < n_rows; i++) {
		inv = inv || (m.getColumn(n_rows + i));
	}
	return inv;
}

Matrix Matrix::operator~() {
	Matrix m(n_cols, n_rows);
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			m.data[j][i] = data[i][j];
		}
	}
	return m;
}




bool Matrix::isOrthogonal() {
	Matrix m = ~(*this);
	Matrix k(m * (*this));
	return identityMatrix(k.n_rows) == k;
}

bool Matrix::isSymmetric() {
	return *this == ~(*this);
}

bool Matrix::isLinearlyIndependent() {
	Basis b(*this);
	return b.dimensions() == n_cols;
}

bool Matrix::getQR(Matrix& q, Matrix& r) {
	if (!isLinearlyIndependent()) {
		return false;
	}
	q = Matrix(n_rows, 0);
	r = Matrix(n_cols, n_cols);
	r.set(!getColumn(0), 0, 0);
	q = q || (getColumn(0) / r.get(0, 0));
	for (unsigned int c = 1; c < n_cols; c++) {
		Vector vPer = getColumn(c);
		for (unsigned int rr = 0; rr <= c; rr++) {
			if (rr != c) {
				r.set(q.getColumn(rr) * getColumn(c), rr, c);
				vPer -= q.getColumn(rr) * r.get(rr, c);
			}
			else {
				r.set(!vPer, rr, c);
				q = q || (vPer / r.get(rr, c));
			}
		}
	}
	return true;
}

Vector Matrix::getLSS(Vector b) {
	//if (rank() == n_cols) {
		return (!((~(*this)) * (*this))) * (~(*this)) * b;
	//}
	return b;
}

Matrix Matrix::getOrthoProj() {
	Matrix AT = ~(*this);
	Matrix A = *this;
	return A * (!(AT * A)) * AT;
}

Complex Matrix::determinant() {
	if (n_rows == n_cols) {
		Matrix m = *this;
		//m.print();
		Complex determinant = 1;
		int offset = 0;
		for (int i = 0; i < m.n_rows; i++) {
			while ((i + offset) < n_cols && m.data[i][i + offset] == 0) { offset++; }
			if (m.data[i][i + offset] != 0) {
				determinant *= m.data[i][i + offset];
				m.multRow(Complex(1) / m.data[i][i + offset], i);
				for (int b = 0; b < n_rows; b++) {
					if (b != i) {
						m.addMultRowTo(-m.data[b][i + offset] / m.data[i][i + offset], i, b);
					}
				}
			}
			else {
				break;
			}
		}
		for (int i = 0; i < n_rows; i++) {
			determinant *= m.data[i][i];
		}
		return determinant;
	}
	return 0;
}

bool Matrix::isPerpendicular(Vector b) {
	Matrix k = ~(*this);
	return (k * b).isZeroVector();
}

Basis Matrix::getImage() {
	return Basis(*this);
}

Basis Matrix::getNullImage() {
	Basis b(*this);
	return !b;
}

Complex Matrix::getVolume() {
	Matrix m(*this);
	Matrix mt = ~m;

	return ((mt * m).determinant()).sqrt();
}

bool Matrix::swapRows(unsigned int a, unsigned int b) {
	if (a < 0 || b < 0 || a >= n_rows || b >= n_rows) {
		return false;
	}
	Vector a_row = getRow(a);
	data[a] = getRow(b);
	data[b] = a_row;
	return true;
}

Matrix Matrix::getCofactor(unsigned int r, unsigned int c) {
	Matrix temp_m(0, n_cols);
	for (int i = 0; i < n_rows; i++) {
		if (i != r) {
			temp_m = temp_m && getRow(i);
		}
	}

	Matrix result(n_rows - 1, 0);
	for (int i = 0; i < n_cols; i++) {
		if (i != c) {
			result = result || temp_m.getColumn(i);
		}
	}

	return result;
}

Matrix Matrix::getAdjoint() {
	Matrix m(n_rows, n_rows);
	for (int r = 0; r < n_rows; r++) {
		for (int c = 0; c < n_cols; c++) {
			m.set(getCofactor(r, c).determinant() * std::pow(-1, r + c), r, c);
		}
	}
	return ~m;
}



Vector Matrix::getRow(int i) const {
	if (i < n_rows) {
		return Vector(data[i]);
	}
	return Vector();
}

Vector Matrix::getColumn(int b) const {
	if (b < n_cols) {
		Vector v;
		for (int i = 0; i < n_rows; i++) {
			v.push(data[i].get(b));
		}
		return v;
	}
	return Vector();
}

void Matrix::setRow(Vector v, int i) {
	while (i >= n_rows) {
		*this = *this && Vector(n_cols);
	}
	data[i] = v;
}

void Matrix::setColumn(Vector v, int i) {
	while (i >= n_cols) {
		*this = *this || Vector(n_rows);
	}
	for (int b = 0; b < v.size(); b++) {
		data[b][i] = v[b];
	}
}





void Matrix::print() {
	cout << endl;
	for (int i = 0; i < n_rows; i++) {
		for (int b = 0; b < n_cols; b++) {
			cout << data[i][b].toString() << " ";
		}
		cout << endl;
	}
	if (n_rows == 0 || n_cols == 0) {
		cout << "< Empty Matrix: " << n_rows << " rows by " << n_cols << " cols >";
	}
	cout << endl;

}


unsigned int Matrix::getLowestNonZeroRow() {
	int i = n_rows - 1;
	while (getRow(i).isZeroVector()) {
		i--;
	}
	return i;
}

bool Matrix::isIdentity() {
	bool ans = true;
	for (int r = 0; r < n_rows; r++) {
		for (int c = 0; c < n_cols; c++) {
			if (r != c) {
				ans = ans && get(r, c) == 0;
			}
			else {
				ans = ans && get(r, c) == 1;
			}
		}
	}
	return ans;
}

Complex Matrix::trace() {
	Complex c;
	for (int i = 0; i < n_rows; i++) {
		c += data[i][i];
	}
	return c;
}

double factorial(int x) {
	if (x < 0) { return -1; }
	if (x <= 1) { return 1; }
	return x * factorial(x - 1);
}

Matrix Matrix::exp() {

	Matrix I(n_rows, n_cols);
	Matrix sol(n_rows, n_cols);

	for (int i = 0; i < n_rows; i++) {
		I.set(Complex(1), i, i);
	}


	for (int i = 0; i < 25; i++) {

		Complex c = Complex(1 / (double)factorial(i));

		Matrix step = I * pow(i);

		step = step * c;

		sol = sol + step;

	}

	return sol;

}

Matrix Matrix::pow(int p) {

	Matrix I(n_rows, n_cols);
	Matrix m(*this);
	for (int i = 0; i < n_rows; i++) {
		I.set(Complex(1), i, i);
	}

	for (int i = 0; i < p; i++) {
		I = I * m;
	}

	return I;

}

Complex Matrix::rightDiag(int start) {
	int c = start;
	Complex result(1, 0);
	for (int r = 0; r < n_rows; r++) {
		result *= data[r][c];
		c = (c + 1) % n_cols;
	}
	return result;
}
Complex Matrix::leftDiag(int start) {
	int c = start;
	Complex result(1, 0);
	for (int r = 0; r < n_rows; r++) {
		result *= data[r][c];
		c = (c - 1) < 0 ? n_cols - 1 : c - 1;
	}
	return result;
}

/*Complex lambdaDet(Matrix m, Complex lambda) {
	for (int i = 0; i < m.numCols(); i++) {
		m.set(((m.get(i, i) - lambda) * 1000000).round() / 1000000, i, i);
	}
	return m.determinant();
}

Complex zoomLambda(Matrix m, Complex lambmin, Complex detmin, Complex lambmax, Complex detmax, int steps) {
	if (detmax < detmin) {
		Complex dm = detmax;
		Complex lm = lambmax;
		lambmax = lambmin;
		detmax = detmin;
		lambmin = lm;
		detmin = dm;
	}
	
	if (steps <= 0) {
		return (lambmin + lambmax) / 2;
	}

	Complex lambmid = (lambmin + lambmax) / 2;
	Complex detmid = lambdaDet(m, lambmid);

	if (detmid < 0) {
		return zoomLambda(m, lambmid, detmid, lambmax, detmax, steps - 1);
	}
	if (detmid > 0) {
		return zoomLambda(m, lambmin, detmin, lambmid, detmid, steps - 1);
	}
	return lambmid;
}

Complex zoomLambdaDeriv(Matrix m, Complex lambmin, Complex detmin, Complex lambmax, Complex detmax, int steps) {

	if (steps <= 0) { return (lambmax + lambmin) / 2; }

	Complex stepsize = (lambmax - lambmin) / 100;
	Complex detip = lambdaDet(m, lambmin - stepsize);
	Complex minDeriv = ((lambdaDet(m, lambmin) - detip) / stepsize);
	Complex minLambda = lambmin - stepsize;

	//cout << lambmin << " " << lambmax << " " << stepsize << endl;

	for (Complex i = lambmin; i <= lambmax; i += stepsize) {
		Complex deti = lambdaDet(m, i);
		Complex deriv = (deti - detip) / stepsize;
		if (deriv.mag() < minDeriv.mag()) {
			minDeriv = deriv.mag();
			minLambda = i;
		}
		detip = deti;
	}

	Complex lm = minLambda - stepsize;
	Complex dm = lambdaDet(m, lm);
	Complex lx = minLambda + stepsize;
	Complex dx = lambdaDet(m, lx);

	return zoomLambdaDeriv(m, lm, dm, lx, dx, steps - 1);

}

bool crosses(Complex a, Complex b) {
	return (a * b) < 0;
}*/


vector<Complex> Matrix::getEigenValues() {
	

	Matrix m(*this);
	vector<Complex> v;
	if (n_rows == 2 && n_cols == 2) {
		Complex a(1);
		Complex b = m.trace();
		Complex c = m.determinant();

		Complex root_1 = (-b + (b * b - a * c * 4).sqrt()) / (a * 2);
		Complex root_2 = (-b - (b * b - a * c * 4).sqrt()) / (a * 2);

		v.push_back(root_2);
		v.push_back(root_1);
	}


	if (n_rows == 3 && n_cols == 3) {
		Complex a(data[0][0]);
		Complex b(data[0][1]);
		Complex c(data[0][2]);
		Complex d(data[1][0]);
		Complex e(data[1][1]);
		Complex f(data[1][2]);
		Complex g(data[2][0]);
		Complex h(data[2][1]);
		Complex i(data[2][2]);

		Complex aa(-1);
		Complex bb(a + e + i);
		Complex cc(f * h + b * d + c * g - e * i - e * a - a * i);
		Complex dd(m.determinant());

		Complex p(-bb / (aa * 3));
		Complex q(p.pow(3) + (bb * cc - aa * dd * 3) / (a * a * 6));
		Complex r(cc / (aa * 3));

		Complex root_1((q + (q * q + (r - p * p).pow(3)).sqrt()).pow(1.0 / 3.0) + (q - (q * q + (r - p * p).pow(3)).sqrt()).pow(1.0 / 3.0) + p);

		Complex aaa(aa);
		Complex bbb(bb + aa * root_1);
		Complex ccc(cc + bb * root_1 + aa * root_1 * root_1);

		Complex root_2 = (-bbb + (bbb * bbb - aaa * ccc * 4).sqrt()) / (aaa * 2);
		Complex root_3 = (-bbb - (bbb * bbb - aaa * ccc * 4).sqrt()) / (aaa * 2);

		v.push_back(root_2.softRound());

		if (root_3 != root_2) {
			v.push_back(root_3.softRound());
		}
		if (root_1 != root_2 && root_1 != root_3) {
			v.push_back(root_1.softRound());
		}
		
	}

	return v;
}




vector<EigenVector> Matrix::getEigenVectors() {

	vector<Complex> eigenvals;
	eigenvals.push_back(Complex(-1));
	eigenvals.push_back(Complex(8));
	vector<EigenVector> results;

	for (int i = 0; i < eigenvals.size(); i++) {
		Complex eigenval = eigenvals[i];
		Matrix m(*this);
		eigenval.print();
		for (int i = 0; i < m.numCols(); i++) {
			m.set(((m.get(i, i) - eigenval) * 1000000).round() / 1000000, i, i);
		}
		m.print();
		for (int i = 0; i < m.numRows(); i++) {
			m.set(-m.get(i, n_rows - 1), i, n_rows - 1);
		}
		m.print();
		m.toRREF();
		m.print();

		Vector NormalAnswer(n_rows);

		for (int a = 0; a < n_rows - 1; a++) {
			for (int b = 0; b < n_rows - 1; b++) {
				if (m.get(a,b) != 0) {
					NormalAnswer[b] = m.get(a, n_cols - 1);
				}
			}
		}

		NormalAnswer[n_cols - 1] = 1;
		NormalAnswer.wholeNumberize();

		results.push_back(EigenVector(eigenval, NormalAnswer));

	}
		

	/*
	
	for (int i = 0; i < eigenvals.size(); i++) {
		Complex eigenval = eigenvals[i];
		Matrix m(*this);
		for (int i = 0; i < m.numCols(); i++) {
			m.set(((m.get(i, i) - eigenval) * 1000000).round()/1000000, i, i);
		}

		m.toRREF();
		m.print();
		vector<int> deletedCols;

		bool foundZeroCol = false;
		for (int i = 0; i < m.numCols(); i++) {
			if (m.getColumn(i).isZeroVector()) {
				deletedCols.push_back(i);
				cout << i << endl;
				m.delColumn(i);
				m.delRow(m.numRows() - 1);
				Vector v(numCols());
				v[i] = 1;
				results.push_back(EigenVector(eigenval, v));
				i--;
			}
		}

		if (m.isIdentity()){
			continue;
		}
		
		int n = m.numCols();
		Vector ans(n);
		Vector v = Vector(n);
		v[n - 1] = 1;
		while ((v - m.getRow(n - 2)).isZeroVector()) {
			m.delRow(n - 2);
			m.delColumn(n - 1);
			n--;
			v = Vector(n);
			v[n - 1] = 1;
		}

		ans[n - 1] = 1;

		for (int i = n - 2; i >= 0; i--) {
			ans[i] = -m.get(i, n - 1);
			if ((ans[i]).mag() < 0.0000001) { ans[i] = 0; }
		}
		ans.wholeNumberize();
		for (int i = 0; i < deletedCols.size(); i++) {
			ans.insert(0, deletedCols[i]);
		}
		results.push_back(EigenVector(eigenval, ans));

	}*/

	return results;
}



Matrix identityMatrix(unsigned int n) {
	Matrix m(n, n);
	for (int i = 0; i < n; i++) {
		m.set(1, i, i);
	}
	return m;
}

Complex CramersRule(Matrix m, Vector sol, unsigned int var) {
	Matrix mk = m;
	mk.setColumn(sol, var);
	return mk.determinant() / m.determinant();
}

Vector CramersRule(Matrix m, Vector sol) {
	Vector v(sol.size());
	for (int i = 0; i < sol.size(); i++) {
		v[i] = CramersRule(m, sol, i);
	}
	return v;
}

Vector findEigenvector(Matrix m) {
	Complex mag = !m.getRow(0);

	for (int i = 0; i < 2; i++) {
			m.set(m.get(0, i) / mag, 0, i);
	}
	m.set(m.get(0 ,0) - 1, 0, 0);

	Vector v(2);
	v[0] = 1;
	v[1] = -m.get(0, 0) / m.get(0, 1);

	return v;
}