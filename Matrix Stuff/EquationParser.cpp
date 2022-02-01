#include "EquationParser.h"



vector<string> split(const string& str, char delim = ' ')
{
	vector<string> cont;
	std::stringstream ss(str);
	std::string token;
	while (std::getline(ss, token, delim)) {
		cont.push_back(token);
	}
	return cont;
}

vector<string> splitSigns(const string& str)
{
	vector<string> cont;
	int index = 0;
	int parenDeep = 0;
	cont.push_back("");
	for (int i = 0; i < str.length(); i++) {
		if (str[i] == '(') { parenDeep++; }
		else if (str[i] == ')') { parenDeep--; }
		
		if (parenDeep == 0) {
			if (str[i] == '+' && cont[index].length() > 0) {
				index++;
				cont.push_back("");
			}
			else if (str[i] == '+') {
				continue;
			}
			else if (str[i] == '-') {
				if (cont[index].length() == 0) {
					cont[index] += str[i];
				}
				else {
					index++;
					cont.push_back("-");
				}
			}
			else {
				cont[index] += str[i];
			}
		}
		else {
cont[index] += str[i];
		}
	}
	return cont;
}

void removeSpaces(string& s) {
	string newstr = "";
	for (int i = 0; i < s.length(); i++) {
		if (s[i] != ' ') { newstr += s[i]; }
	}
	s = newstr;
}

vector<double> parseEquation(string s) { return vector<double>(); };

struct Term {
	Term(string s) {
		vector<string> temp = split(s, '*');
		mult = stod(temp[0]);
		var = temp[1][0];
	}
	double mult;
	char var;
};

vector<double> simpleEval(string s, string varNames) {

	//cout << s << endl;

	removeSpaces(s);

	vector<string> split_ = split(s, '=');

	vector<string> leftSideTemp = splitSigns(split_[0]);

	if (leftSideTemp.size() > varNames.length()) {
		cout << "Too Many Variables in Equation!" << endl;
		return vector<double>();
	}

	vector<Term> leftSide;

	for (int i = 0; i < leftSideTemp.size(); i++) {
		leftSide.push_back(Term(leftSideTemp[i]));
	}

	vector<double> solution;

	for (int i = 0; i < varNames.length(); i++) {
		int index = -1;
		for (int b = 0; b < leftSide.size(); b++) {
			if (varNames[i] == leftSide[b].var) {
				index = b;
			}
		}
		if (index != -1) {
			solution.push_back(leftSide[index].mult);
		}
		else {
			solution.push_back(0);
		}
	}

	if (split_.size() > 1) {

		solution.push_back(stod(split_[1]));

	}

	else {

		solution.push_back(0);

	}


	return solution;


}