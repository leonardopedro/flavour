#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>
using namespace std;
using namespace Eigen;

class points{
public:
	points(double p, double m, double e):prediction(p),measure(m),error(e){}
double prediction, measure, error;	
};

double chi2(double * v){
	double buu = v[0], auu = v[1];
	double add = v[2], bdd = v[3];
	double eu =v[4], gu = v[5];
	double ed = v[6], gd = v[7];
complex<double> bu = buu*exp(complex<double>(0,gu*M_PI_2)), au = auu*exp(complex<double>(0,eu*M_PI_2));
complex<double> ad = add*exp(complex<double>(0,ed*M_PI_2)), bd = bdd*exp(complex<double>(0,gd*M_PI_2)); 
complex<double> cd = bd - ad; 
complex<double> cu =bu + au;

//Matrix3cd X = Matrix3cd::Random(3,3);
Matrix3cd mu,md;
mu<<1,1,1.0+au+bu,1,1,1.0+bu,1.0+bu+au,1.0+bu,1.0+cu;
md<<1,1,1.0+ad+bd,1,1,1.0+bd,1.0+bd+ad,1.0+bd,1.0+cd;
	mu=mu*v[8];
	md=md*v[9];

Matrix3cd Hu = mu*mu.adjoint(), Hd = md*md.adjoint();
SelfAdjointEigenSolver<Matrix3cd> usolver(Hu), dsolver(Hd);
Vector3d Du=usolver.eigenvalues();//.asDiagonal();
Vector3d Dd=dsolver.eigenvalues();//.asDiagonal();
Matrix3cd VLd=dsolver.eigenvectors().adjoint();
Matrix3cd VLu=usolver.eigenvectors().adjoint();
Matrix3cd Vckm=VLu*VLd.adjoint();
double lambda=sqrt(norm(Vckm(0,1))/(norm(Vckm(0,0))+norm(Vckm(0,1))));
double A=abs(Vckm(1,2))/abs(Vckm(0,1))/lambda;
complex<double> rhoeta=-Vckm(0,0)*conj(Vckm(0,2))/Vckm(1,0)/conj(Vckm(1,2));
double rho=real(rhoeta), eta=imag(rhoeta);

vector<points> comp;
comp.push_back(points(lambda, 0.22535, 0.00065));
comp.push_back(points(A, 0.811, 0.017));
comp.push_back(points(rho, 0.131, 0.02));
comp.push_back(points(eta, 0.345, 0.014));
comp.push_back(points(eta, 0.345, 0.014));
comp.push_back(points(Du[0], 1.27e-3, 0.46e-3));
comp.push_back(points(Dd[0], 2.9e-3, 1.22e-3));
comp.push_back(points(Dd[1], 55e-3, 16e-3));
comp.push_back(points(Du[1], 0.619, 0.084));
comp.push_back(points(Dd[2], 2.89, 0.09));
comp.push_back(points(Du[2], 171.7, 3.0));

double chi=0;

for(uint i=0;i<comp.size();i++){
	chi+=pow((comp[i].prediction-comp[i].measure)/comp[i].error,2);
	}
return chi/comp.size();
}

int main(){
double buu = 0.0162, auu = 0.0009;
double add = 0.018, bdd = 0.09;
//double x = 1.0/3;
double eu = -1.8, gu = -1.0/2;
double ed = -1.0/2, gd = -2;
complex<double> bu = buu*exp(complex<double>(0,gu*M_PI_2)), au = auu*exp(complex<double>(0,eu*M_PI_2));
complex<double> ad = add*exp(complex<double>(0,ed*M_PI_2)), bd = bdd*exp(complex<double>(0,gd*M_PI_2)); 
complex<double> cd = exp(complex<double>(0,0.09/1.8)); 

//Matrix3cd X = Matrix3cd::Random(3,3);
Matrix3cd mu,md;
mu<<1,1,1.0+au+bu,1,1,1.0+bu,1.0+bu+au,1.0+bu,1.0+bu;
//md<<1,1,1.0+ad+bd,1,1,1.0+bd,1.0+bd+ad,1.0+bd,1.0+bd;
md<<1,1,1.0+ad+bd,1,1,1.0+bd,cd*(1.0+bd+ad),cd*(1.0+bd),cd*(1.0+bd);
mu=mu*57.4822;
md=md*1.0147;

Matrix3cd Hu = mu*mu.adjoint(), Hd = md*md.adjoint();
SelfAdjointEigenSolver<Matrix3cd> usolver(Hu), dsolver(Hd);
Vector3d Du=usolver.eigenvalues();//.asDiagonal();
Vector3d Dd=dsolver.eigenvalues();//.asDiagonal();
Matrix3cd VLd=dsolver.eigenvectors().adjoint();
Matrix3cd VLu=usolver.eigenvectors().adjoint();
Matrix3cd Vckm=VLu*VLd.adjoint();
double lambda=sqrt(norm(Vckm(0,1))/(norm(Vckm(0,0))+norm(Vckm(0,1))));
double A=abs(Vckm(1,2))/abs(Vckm(0,1))/lambda;
complex<double> rhoeta=-Vckm(0,0)*conj(Vckm(0,2))/Vckm(1,0)/conj(Vckm(1,2));
double rho=real(rhoeta), eta=imag(rhoeta);
double fu=173.5/sqrt(Du[2]), fd=2.89/sqrt(Dd[2]);

vector<points> comp;
comp.push_back(points(lambda, 0.22535, 0.00065));
comp.push_back(points(A, 0.811, 0.017));
comp.push_back(points(rho, 0.131, 0.02));
comp.push_back(points(eta, 0.345, 0.014));
comp.push_back(points(sqrt(Du[0]), 1.27e-3, 0.46e-3));
comp.push_back(points(sqrt(Dd[0]), 2.9e-3, 1.22e-3));
comp.push_back(points(sqrt(Dd[1]), 55e-3, 16e-3));
comp.push_back(points(sqrt(Du[1]), 0.619, 0.084));
comp.push_back(points(sqrt(Dd[2]), 2.89, 0.09));
comp.push_back(points(sqrt(Du[2]), 171.7, 3.0));

double chi=0;

for(uint i=0;i<comp.size();i++){
	double aa=(comp[i].prediction-comp[i].measure)/comp[i].error;
	cout<<i<<" "<<aa<<endl;
	chi+=pow((comp[i].prediction-comp[i].measure)/comp[i].error,2);
	}
chi/=comp.size();

cout<<lambda<<" "<<A<<" "<<rho<<" "<<eta<<endl;
cout<<sqrt(Du[2])<<" "<<sqrt(Du[1])<<" "<<sqrt(Du[0])<<endl;
cout<<sqrt(Dd[2])<<" "<<sqrt(Dd[1])<<" "<<sqrt(Dd[0])<<endl;
cout<<chi<<endl;
/*
cout << "The eigenvalues of mu are:" << endl << es.eigenvalues() << endl;
cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
double lambda = es.eigenvalues()[0];
cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
Vector3cd v = es.eigenvectors().col(0);
cout << "If v is the corresponding eigenvector, then lambda * v = " << endl << lambda * v << endl;
cout << "... and A * v = " << endl << A * v << endl << endl;
Matrix3d D = es.eigenvalues().asDiagonal();
Matrix3cd V = es.eigenvectors();
cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
*/
return 0;
}
