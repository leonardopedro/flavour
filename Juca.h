#ifndef JUCA_H
#define JUCA_H

#include "widthcalc.h"

#include "TH2F.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include <iostream>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

class Juca: public Model{
public:

Juca():Mu("Mu"), Md("Md"), Mc("Mc"), Ms("Ms"), Mt("Mt"), Mb("Mb"), lambda("lambda"), A("A"),rho("rho"),eta("eta"){
    
	add("",lambda,new gauss2obs(0.22535,0.00065));
	add("",A,new gauss2obs(0.811, 0.017));
	add("",rho,new gauss2obs(0.131, 0.02));
	add("",eta,new gauss2obs(0.345, 0.014));
	add("",Mu,new gauss2obs(1.27e-3, 0.46e-3));
	add("",Md,new gauss2obs(2.9e-3, 1.22e-3));
	add("",Ms,new gauss2obs(55e-3, 16e-3));
	add("",Mc,new gauss2obs(0.619, 0.084));
	add("",Mb,new gauss2obs(2.89, 0.09));
	add("",Mt,new gauss2obs(171.7, 3.0));
}

~Juca(){}



void add(const char * s, ex pred, observable * ob){
	ex p=collect_common_factors(expand(pred.subs(replacements)));
	push_back(prediction(ob,p));
	}
	
parameters generateparameters() const{
	parameters p;
	
	//buu
	p.push_back(freeparameter(-5,0,r));
	//auu
	p.push_back(freeparameter(-5,0,r));
	//add
	p.push_back(freeparameter(-5,0,r));
	//bdd
	p.push_back(freeparameter(-5,0,r));
	//eu
	p.push_back(freeparameter(-2,2,r));
	//gu
	p.push_back(freeparameter(-2,2,r));
	//ed
	p.push_back(freeparameter(-2,2,r));
	//gd
	p.push_back(freeparameter(-2,2,r));
	//nu
	p.push_back(freeparameter(-2,3,r));
	//nd
	p.push_back(freeparameter(-2,3,r));
	//cuu
	//p.push_back(freeparameter(-5,0,r));
	//cdd
	p.push_back(freeparameter(-5,0,r));
	//hu
	//p.push_back(freeparameter(-2,2,r));
	//hd
	p.push_back(freeparameter(-2,2,r));
	
	return p;
}


lst getlist(const parameters & p) const{
	
	
	double buu = pow(10.0,p[0].value), auu = pow(10.0,p[1].value);
	double add = pow(10.0,p[2].value), bdd = pow(10.0,p[3].value);
	double eu =p[4].value, gu = p[5].value;
	double ed = p[6].value, gd = p[7].value;
	//double cuu = pow(10.0,p[10].value); 
	double cdd = pow(10.0,p[10].value);
	//double hu =p[12].value;
	double hd= p[11].value;
	
complex<double> bu = buu*exp(complex<double>(0,gu*M_PI_2)), au = auu*exp(complex<double>(0,eu*M_PI_2));
complex<double> ad = add*exp(complex<double>(0,ed*M_PI_2)), bd = bdd*exp(complex<double>(0,gd*M_PI_2)); 
//complex<double> cu = cuu*exp(complex<double>(0,hu*M_PI_2));
complex<double> cu = bu;
complex<double> cd = cdd*exp(complex<double>(0,hd*M_PI_2));
//complex<double> cd = bd;

//Matrix3cd X = Matrix3cd::Random(3,3);
Matrix3cd mu,md;
mu<<1,1.0,1.0+au+bu,1.0,1.0,1.0+bu,1.0+bu+au,1.0+bu,1.0+cu;
md<<1,1.0,1.0+ad+bd,1.0,1.0,1.0+bd,1.0+bd+ad,1.0+bd,1.0+cd;
	mu=mu*pow(10.0,p[8].value);
	md=md*pow(10.0,p[9].value);

Matrix3cd Hu = mu*mu.adjoint(), Hd = md*md.adjoint();
SelfAdjointEigenSolver<Matrix3cd> usolver(Hu), dsolver(Hd);
Vector3d Du=usolver.eigenvalues();//.asDiagonal();
Vector3d Dd=dsolver.eigenvalues();//.asDiagonal();
Matrix3cd VLd=dsolver.eigenvectors().adjoint();
Matrix3cd VLu=usolver.eigenvectors().adjoint();
Matrix3cd Vckm=VLu*VLd.adjoint();
double lambda0=sqrt(norm(Vckm(0,1))/(norm(Vckm(0,0))+norm(Vckm(0,1))));
double A0=abs(Vckm(1,2))/abs(Vckm(0,1))/lambda0;
complex<double> rhoeta=-Vckm(0,0)*conj(Vckm(0,2))/Vckm(1,0)/conj(Vckm(1,2));
double rho0=real(rhoeta), eta0=imag(rhoeta);
	
	return lst(lambda==lambda0,A==A0,rho==rho0,eta==eta0,Mu==sqrt(abs(Du[0])),Md==sqrt(abs(Dd[0])),Mc==sqrt(abs(Du[1])),Ms==sqrt(abs(Dd[1])),Mt==sqrt(abs(Du[2])),Mb==sqrt(abs(Dd[2])));
}

const possymbol Mu, Md, Mc, Ms, Mt, Mb, lambda, A;
const realsymbol rho, eta;

lst replacements;
};



#endif
