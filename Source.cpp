#include <iostream> 
#include <math.h>

using namespace std;

const double g = 9.8; // ускорение свободного падения
const double PI = 3.14;

double func(double x, double u, double S, double R) {			// Исходная функция du/dx=func()
	return -0.6 * S * pow(2 * g, 0.5) * pow(PI ,-1) * pow(u, -0.5) * pow(2 * R - u, -1);
}

double k1(double x, double u, double S, double R, double h) {
	return func(x, u, S, R);
}

double k2(double x, double u, double S, double R, double h, double k1) {
	return func(x + h / 3, u + k1 * h / 3, S, R);
}

double k3(double x, double u, double S, double R, double h, double k1, double k2) {
	return func(x + 2 * h / 3, u - h * k1 / 3 + h * k2, S, R);
}

double k4(double x, double u, double S, double R, double h, double k1, double k2, double k3) {
	return func(x + h, u + h * k1 - h * k2 + h * k3, S, R);
}

double solve(double x0, double u0, double El, double S, double R) {
	double u = u0;
	double x = x0;
	cout << "x = " << x0 << "\t\t" << "u = " << u0 << endl;

	do {
		double h = 0.02;
		double uh;
		double uhp;

		do {
			h /= 2;
			// решение c шагом h
			double k1_ = k1(x, u, S, R, h);
			double k2_ = k2(x, u, S, R, h, k1_);
			double k3_ = k3(x, u, S, R, h, k1_, k2_);
			double k4_ = k4(x, u, S, R, h, k1_, k2_, k3_);

			uh = u + h / 8 * (k1_ + 3 * k2_ + 3 * k3_ + k4_);

			// решение с шагом h/2
			k1_ = k1(x, u, S, R, h/2);
			k2_ = k2(x, u, S, R, h/2, k1_);
			k3_ = k3(x, u, S, R, h/2, k1_, k2_);
			k4_ = k4(x, u, S, R, h/2, k1_, k2_, k3_);

			uhp = u + h / 16 * (k1_ + 3 * k2_ + 3 * k3_ + k4_);

			k1_ = k1(x + h / 2, uhp, S, R, h / 2);
			k2_ = k2(x + h / 2, uhp, S, R, h / 2, k1_);
			k3_ = k3(x + h / 2, uhp, S, R, h / 2, k1_, k2_);
			k4_ = k4(x + h / 2, uhp, S, R, h / 2, k1_, k2_, k3_);

			uhp = uhp + h / 16 * (k1_ + 3 * k2_ + 3 * k3_ + k4_);
			
		} while (abs((uh - uhp) / (1 - (double)1 / 16)) > El);

		u = uh;
		x += h;
		
		if (uh >= 0 && uhp >= 0)
			cout << "x = " << x << "\t\t" << "u = " << uh << "\t\tЛокальная погрешность: " << abs((uh - uhp) / (1 - (double)1 / 16)) << endl;



	} while (u >= 0);

	return x;
}

int main() {
	setlocale(LC_ALL, "rus");
	double El = 10e-10; // локальная погрешность
	double x0 = 0;		// начальное время
	double u0 = 0.1;     // u(0) = u0 начальная высота
	double S = 0.1;		// площадь отверстия для слива
	double R = 0.2;		// радиус

	double result = solve(x0, u0, El, S, R);

	cout << "Result: " << result;

	return 0;
}