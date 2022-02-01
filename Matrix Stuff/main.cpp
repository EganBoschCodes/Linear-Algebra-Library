#include <cmath>
#include "Matrix.h"
#include "EquationParser.h"
#include "Quantity.h"
#include "Vector.h"
#include "Basis.h"
#include "EquationSolver.h"
#include "Complex.h"

using namespace std;




int floor_(double d) {
	return (int)d;
}


double yprime(double x, double y) {
	return x * x - x * y;
}



int main() {

	/*Complex d = -2;
	d.print();

	cout << d.mag() << ", " << d.angle() << endl;

	Complex c = d.sqrt();
	c.print();

	cout << c.mag() << ", " << c.angle() << endl;*/

	Vector u1("2 -3 4 -5 2");

	Vector u3("3 -2 7 -9 1");

	Vector u5("-1 1 2 1 -3");

	Vector u7("1 0 -2 3 -2");

	Vector u8("2 -1 1 -9 7");

	Matrix a(5, 0);

	a = a || u1 || u3 || u5 || u7 || u8;

	a.print();

	a.toRREF();

	a.print();

	Vector test = (u1 * 4) - (u3 * 2);
	test.print();


	/*Matrix a("3 2 4,2 0 2,4 2 3");

	vector<EigenVector> v = a.getEigenVectors();
	for (int i = 0; i < v.size(); i++) {
		v[i].print();
	}
	Vector v2 = eigens[eigens.size()-1].vector();
	Vector v3 = m * v2;
	v2.print();
	v3.print();*/


}




//Equation Solver

/*int main() {
	while (true) {
		int num_vars;
		string var_string;
		cout << "What variables will you be working with (written together, no commas or spaces, like \"xyz\"): ";
		getline(cin, var_string);
		num_vars = var_string.length();


		Matrix m(num_vars, num_vars + 1);
		string equation;
		string equations[3] = { "y = 2*x-1", "z = -0.5*x + 1.5", "x = 0" };
		for (int i = 0; i < num_vars; i++) {
			cout << "Equation " << (i + 1) << ": " << equations[i];
			m.setRow(Equation(equations[i]).eval(var_string), i, 0);
			cout << endl;
		}
		cout << "Matrix before solution: " << endl;
		m.print();


		m.toRREF();
		cout << "Matrix after solution: " << endl;
		m.print();


		cout << "Solution: " << endl;
		for (int i = 0; i < num_vars; i++) {
			cout << var_string[i] << " = " << m.get(i, num_vars) << endl;
		}
		cout << endl << endl;
	}

	return 0;
}*/
