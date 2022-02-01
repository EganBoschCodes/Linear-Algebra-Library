#include "Quantity.h"
#include "EquationParser.h"
#include <string>

using namespace std;








Quantity::Quantity(std::string q) {


	removeSpaces(q);
	while (q[0] == '(' && q[q.length() - 1] == ')') { 
		q = q.substr(1, q.length() - 2); 
	}
	//cout << q << endl;
	k = 0;
	vector<string> split_sign = splitSigns(q);
	if (split_sign.size() == 1) {
		size_t parenpos = q.find('(');
		if (parenpos == string::npos) {
			vector<string> data = split(split_sign[0], '*');
			if (data.size() == 1) {
				if (data[0].find("/") != string::npos) {
					data = split(data[0], '/');
					if (isalpha(data[0][0])) {
						add(1, data[0][0]);
						mult(1 / stod(data[1]));
						return;
					}
					else {
						add(stod(data[0]) / stod(data[1]));
					}
				}
				else {
					if (isalpha(data[0][0])) {
						add(1, data[0][0]);
					}
					else if (isalpha(data[0][1]) && data[0][0] == '-') {
						add(-1, data[0][1]);
					}
					else {
						add(stod(data[0]));
					}
					return;
				}
			}
			else {
				if (data[0].find("/") != string::npos) {
					vector<string> data2 = split(data[0], '/');
					double mul = stod(data2[0]) / stod(data2[1]);
					add(1, data[1][0]);
					mult(mul);
					return;
				}
				else if (data[1].find("/") != string::npos) {
					vector<string> data2 = split(data[1], '/');
					double mul = stod(data[0]) / stod(data2[1]);
					add(1, data2[0][0]);
					mult(mul);
					return;
				}
				else {
					if (isalpha(data[0][0])) {
						add(1, data[0][0]);
						mult(stod(data[1]));
						return;
					}
					add(1, data[1][0]);
					mult(stod(data[0]));
				}
			}
		}
		else {
			string coeff;
			Quantity parenQuantity = Quantity(q.substr(parenpos, q.length() - parenpos));
			if (q[parenpos - 1] == '*') {
				coeff = q.substr(0, parenpos-1);
			}
			else {
				coeff = q.substr(0, parenpos);
			}
			double d;
			if (coeff.find("/") != string::npos) {
				vector<string> cospl = split(coeff, '/');
				d = stod(cospl[0]) / stod(cospl[1]);
			}
			else {
				d = stod(coeff);
			}
			add(parenQuantity.getMult(d));
		}
	}
	else {
		for (int i = 0; i < split_sign.size(); i++) {
			add(Quantity(split_sign[i]));
		}
	}
	
}

void Quantity::add(ETerm q) {
	ETerm* t = getToAdd(q.var);
	t->val += q.val;
	t->var = q.var;
}

void Quantity::add(Quantity q) {
	k += q.k;
	for (int i = 0; i < q.terms.size(); i++) {
		add(q.terms[i]);
	}
}

void Quantity::mult(double d) {
	k *= d;
	for (int i = 0; i < terms.size(); i++) {
		terms[i].val *= d;
	}
}


Quantity::ETerm* Quantity::getToAdd(char c) {
	for (int i = 0; i < terms.size(); i++) {
		if (terms[i].var == c) {
			return &(terms[i]);
		}
	}
	terms.push_back(ETerm('1', 0));
	return &(terms[terms.size() - 1]);
}

void Quantity::print() {
	cout << toString() << endl;
}

string Quantity::toString() {
	stringstream s;
	for (int i = 0; i < terms.size(); i++) {
		s << terms[i].val << "*" << terms[i].var << (k != 0 || i < terms.size() - 1 ? " + " : "");
	}

	if (k != 0) { s << k; }
	return s.str();
}




Equation::Equation(const std::string eq) {
	vector<std::string> temp = split(eq, '=');

	removeSpaces(temp[0]);
	removeSpaces(temp[1]);

	left = Quantity(Quantity(temp[0]), Quantity(temp[1]).getMult(-1));

	right.add(-left.getConst());
	left.add(-left.getConst());
}

Equation::Equation(Quantity l, Quantity r) {
	left = Quantity(l, r.getMult(-1));
	right.add(-left.getConst());
	left.add(-left.getConst());
}

vector<double> Equation::eval(string vars) {
	return simpleEval(toString(), vars);
}

void Equation::print() {
	cout << toString() << endl;
}

string Equation::toString() {
	stringstream s;
	s << left.toString() << " = " << right.toString();
	return s.str();
}
