#include <iostream>
#include "Eigen/Eigen"
#include <iomanip>

using namespace std;
using namespace Eigen;

double errore_rel(VectorXd x,VectorXd x_esatto)
{
	return (x_esatto-x).norm()/(x_esatto).norm();
}

int main()
{
	Vector2d x_ex {{-1.0e+0},{-1.0e+00}};
	Matrix2d A1 {{5.547001962252291e-01, -3.770900990025203e-02}, {8.320502943378437e-01, -9.992887623566787e-01}};
	Vector2d b1 {{-5.169911863249772e-01}, {1.672384680188350e-01}};
    Matrix2d A2 {{5.547001962252291e-01, -5.540607316466765e-01}, {8.320502943378437e-01, -8.324762492991313e-01}};
	Vector2d b2 {{-6.394645785530173e-04}, {4.259549612877223e-04}};
	Matrix2d A3 {{5.547001962252291e-01, -5.547001955851905e-01}, {8.320502943378437e-01, -8.320502947645361e-01}};
	Vector2d b3 {{-6.400391328043042e-10}, {4.266924591433963e-10}};
	
	Vector2d xl1 = A1.partialPivLu().solve(b1);
	Vector2d xl2 = A2.partialPivLu().solve(b2);
	Vector2d xl3 = A3.partialPivLu().solve(b3);
	
	Vector2d xq1 = A1.householderQr().solve(b1);
	Vector2d xq2 = A2.householderQr().solve(b2);
	Vector2d xq3 = A3.householderQr().solve(b3);
	
	double err_rel1 = errore_rel(xl1,x_ex);
	double err_rel2 = errore_rel(xl2,x_ex);
	double err_rel3 = errore_rel(xl3,x_ex);
	
	double err_relq1 = errore_rel(xq1,x_ex);
	double err_relq2 = errore_rel(xq2,x_ex);
	double err_relq3 = errore_rel(xq3,x_ex);
	
	cout << "1) err. relativo con PA=LU: " << err_rel1 << "    err. relativo con A=QR: " << err_relq1 << endl;
	cout << "2) err. relativo con PA=LU: " << err_rel2 << "    err. relativo con A=QR: " << err_relq2 << endl;
	cout << "3) err. relativo con PA=LU: " << err_rel3 << "    err. relativo con A=QR: " << err_relq3 << endl;
	
	/* vista dei risultati delle tre matrici con entrambe le fattorizzioni
	cout << std::setprecision(16) << std::scientific;
	cout << xl1 << endl << xq1 << endl;
	cout << xl2 << endl << xq2 << endl;
	cout << xl3 << endl << xq3 << endl; */
	
	return 0;
}
