#ifndef EQUATION_SOLVER
#define EQUATION_SOLVER


#include <chrono>
#include <vector>

using namespace std;
using namespace std::chrono;


/*

ZERO FINDER

*/

double zero(double x) {
	return 0;
}

double diff(double (*f1)(double), double (*f2) (double), double val) {
	return f1(val) - f2(val);
}

vector<double> findSolutions(double (*f1)(double), double (*f2) (double), double min, double max, double increment) {
	milliseconds ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());

	vector<double> solutions;

	vector<double> needsRefining;

	//First construct all possible solutions in the frame.
	for (double d = min; d < max; d += increment) {
		double d1 = diff(f1, f2, d);
		if (d1 == 0) {
			solutions.push_back(d);
			continue;
		}
		double d2 = diff(f1, f2, d + increment);
		if ((d1 < 0 && d2 > 0) || (d1 > 0 && d2 < 0)) {
			needsRefining.push_back(d);
		}
	}

	//Then refine the necessary solutions
	for (int i = 0; i < needsRefining.size(); i++) {
		double dx = increment / 1000;
		double x = needsRefining[i];
		//Newton's method
		for (int k = 0; k < 5; k++) {
			double y = diff(f1, f2, x);
			double dy = diff(f1, f2, x + dx) - y;
			double slope = dy / dx;
			x += -y / slope;
			dx /= 10;
		}

		solutions.push_back(x);

	}

	sort(solutions.begin(), solutions.end());

	cout << "Done: Found " << solutions.size() << " solutions between " << min << " and " << max << " with increment " << increment << " in " << (duration_cast<milliseconds>(system_clock::now().time_since_epoch()) - ms).count() << " milliseconds." << endl;

	return solutions;

}

vector<double> findSolutions(double (*f1)(double), double (*f2) (double), double min, double max) {
	return findSolutions(f1, f2, min, max, (max - min) / 10000);
}

vector<double> findSolutions(double (*f1)(double), double (*f2) (double)) {
	return findSolutions(f1, f2, -50, 50);
}

vector<double> findSolutions(double (*f1)(double)) {
	return findSolutions(f1, zero);
}





/*

CALCULUS

*/

double derivative(double (*f)(double), double x) {
	double interval = 1 / 65536;
	return (f(x + interval) - f(x - interval)) / (2 * interval);
}

double integral(double (*f)(double), double x1, double x2, unsigned int intervals) {
	double step = (x2 - x1) / ((double)intervals);
	double answer = 0;
	for (double x = min(x1, x2); x < max(x1, x2); x += abs(step)) {
		answer += f(x) * step;
	}

	return answer;
}

double integral(double (*f)(double), double x1, double x2) {
	return integral(f, x1, x2, (unsigned int)(sqrt(abs(x2 - x1)) * 1000000));
}

double integral(double (*f)(double), double x) {
	return integral(f, 0, x);
}


#endif