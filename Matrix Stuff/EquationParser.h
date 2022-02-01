#ifndef PARSER
#define PARSER

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iterator>

using namespace std;

vector<string> split(const string& str, char delim);

vector<string> splitSigns(const string& str);

void removeSpaces(string& s);


vector<double> parseEquation(string s);

vector<double> simpleEval(string s, string varNames);

#endif