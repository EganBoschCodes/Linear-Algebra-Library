#ifndef EQUATION
#define EQUATION

#include <vector>
#include <string>

using namespace std;


class Quantity {
public:
	Quantity(std::string q);

	Quantity() {
		k = 0;
	}

	Quantity(const Quantity& q) {
		k = q.k;
		terms = q.terms;
	}

	Quantity(const Quantity& q, const Quantity& u) {
		k = q.k;
		terms = q.terms;
		add(u);
	}

	void add(double val, char var) { add(ETerm(var, val)); }
	void add(Quantity q);
	void add(double kk) { k += kk; }
	Quantity getAdd(Quantity q) { return Quantity(*this, q); }

	void mult(double val);
	Quantity getMult(double d) { Quantity q(*this); q.mult(d); return q; }

	double getConst() {
		return k;
	};

	string getVars() {
		string s = "";
		for (int i = 0; i < terms.size(); i++) {
			s += terms[i].var;
		}
		return s;
	};

	void print();
	string toString();

private:
	struct ETerm {
		ETerm(char v, double vl) {
			var = v;
			val = vl;
		}
		char var;
		double val;
	};
	void add(ETerm t);
	ETerm* getToAdd(char c);
	double k;
	std::vector<ETerm> terms;
};

class Equation {
public:
	Equation(std::string eq);
	Equation(Quantity l, Quantity r);

	vector<double> eval(string vars);

	void print();
	string toString();
private:
	Quantity left, right;
};


#endif

