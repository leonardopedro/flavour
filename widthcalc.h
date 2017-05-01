#ifndef WIDTHCALC_H
#define WIDTHCALC_H

#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <set>
#include <ginac/ginac.h>
#include <cmath>
#include <complex>
#include "multivector.h"


//g++ teste.cpp -o teste -lcln -lginac
using namespace std;
using namespace GiNaC;

/**
* @brief this class calculates decay widths of one lepton to 3 leptons  
*/
class widthcalc{

public:

widthcalc(): M2(0,2,2,2,2,2,2,2), M22(0,2,2,2,2,2,2), s2("s2"), s3("s3"), mq1("mq1"), mq2("mq2"), mq3("mq3"), mq4("mq4"){
	
	integral::max_integration_level=100;
	integral::relative_integration_error=1e-3;
	
	genM2();
	genM22();
}

void genM22(){
	cout<<"Generating M22.dat"<<endl;
	
	realsymbol mu("mu"), nu("nu"), alpha("alpha"), beta("beta"), rho("rho"), zeta("zeta");
	realsymbol q1("q1"), q2("q2"), q3("q3"), q4("q4");
	realsymbol h1("h1"), h2("h2"), h3("h3"), h4("h4");
	
	varidx imu(mu,4,0), inu(nu,4,0), irho(rho,4,0);
	varidx ialpha(alpha,4,0), ibeta(beta,4,0), izeta(zeta,4,0);
	
	varidx jmu(mu,4,1), jnu(nu,4,1), jrho(rho,4,1);
	
	ex m2q1=mq1*mq1, m2q2=mq2*mq2, m2q3=mq3*mq3, m2q4=mq4*mq4;
	
	ex q2mu=indexed(q2,imu);
	ex q3mu=indexed(q3,imu);
	ex q4mu=indexed(q4,imu);
	ex q1mu=q2mu+q3mu+q4mu;
	
	ex s4=m2q1+m2q2+m2q3+m2q4-s2-s3;
	
	ex vq1=dirac_slash(q2,4)+dirac_slash(q3,4)+dirac_slash(q4,4)+mq1*dirac_ONE(), vq2=dirac_slash(q2,4)+mq2*dirac_ONE();
	ex vq3=dirac_slash(q3,4)+mq3*dirac_ONE(), vq4=dirac_slash(q4,4)-mq4*dirac_ONE();
	
    scalar_products sp;
    //sp.add(q2, q3, (s4-m2q2-m2q3)/2);
    sp.add(q4, q3, (s2-m2q4-m2q3)/2);
    //sp.add(q2, q4, (s3-m2q2-m2q4)/2);
  
    //sp.add(q2, q2, m2q2);
    sp.add(q3, q3, m2q3);
    sp.add(q4, q4, m2q4);
    
    //sp.add(h1,h1,-1);
    //sp.add(h2,h2,-1);
    //sp.add(h3,h3,-1);
    //sp.add(h4,h4,-1);
    //sp.add(h1,h1,-1);
    
    multivector<ex,3> v(0,2,2,2);
	v[0][0][0]=dirac_gammaL(); v[0][0][1]=dirac_gammaR();
	v[0][1][0]=dirac_gammaR(); v[0][1][1]=dirac_gammaL();
	v[1][0][0]=dirac_gammaL(); v[1][0][1]=dirac_gammaL();
	v[1][1][0]=dirac_gammaR(); v[1][1][1]=dirac_gammaR();
	
	for(uint i=0;i<2;i++)
	  for(uint j=0;j<2;j++){
		v[0][i][j]=-v[0][i][j]*s2/(mq1+mq2);
		v[1][i][j]=(dirac_slash(q3,4)+dirac_slash(q4,4))*v[1][i][j];
	}
	
	//vector<ex> prop(2,0);
	//prop[0]=-s2/(mq1+mq2);
	//prop[1]=indexed(q3,imu.toggle_variance())+indexed(q4,imu.toggle_variance());
	/*
	multivector<ex,2> prop2(0,2,2);
	for(uint i=0;i<2;i++)
	  for(uint j=0;j<2;j++){
				prop2[i][j]=prop[i]*prop[j].subs(mu==nu);
			}
	*/
		
	//ofstream f("M22.dat");
	for(uint i=0;i<2;i++)
	  for(uint j=0;j<2;j++)
	    for(uint k=0;k<2;k++)
			for(uint l=0;l<2;l++)
			for(uint m=0;m<2;m++)
			for(uint n=0;n<2;n++){
				//cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<endl;
				//cout<<dirac_trace(vq3*v[i][m][1]*vq4*v[j][n][0].subs(mu==nu))<<endl;
				
				ex tmp=dirac_trace(vq3*v[i][m][0]*vq4*v[j][n][1])*int(pow(-1.0,double(k+1+l)));
				M22[i][j][k][l][m][n]=tmp.simplify_indexed(sp);
				//cout<<M22[i][j][k][l][m][n]<<endl<<endl;
				//f<<M22[i][j][k][l][m][n]<<endl;
			}
}

void genM2(){
	cout<<"Generating M2.dat"<<endl;
	
	realsymbol mu("mu"), nu("nu"), alpha("alpha"), beta("beta"), rho("rho"), zeta("zeta");
	realsymbol q1("q1"), q2("q2"), q3("q3"), q4("q4");
	realsymbol h1("h1"), h2("h2"), h3("h3"), h4("h4");
		
	varidx imu(mu,4,0), inu(nu,4,0), irho(rho,4,0);
	varidx ialpha(alpha,4,0), ibeta(beta,4,0), izeta(zeta,4,0);
	
	varidx jmu(mu,4,1), jnu(nu,4,1), jrho(rho,4,1);
	
	ex m2q1=mq1*mq1, m2q2=mq2*mq2, m2q3=mq3*mq3, m2q4=mq4*mq4;
	
	ex q2mu=indexed(q2,imu);
	ex q3mu=indexed(q3,imu);
	ex q4mu=indexed(q4,imu);
	ex q1mu=q2mu+q3mu+q4mu;
	
	
	
	ex vq1=dirac_slash(q2,4)+dirac_slash(q3,4)+dirac_slash(q4,4)+mq1*dirac_ONE(), vq2=dirac_slash(q2,4)+mq2*dirac_ONE();
	ex vq3=dirac_slash(q3,4)+mq3*dirac_ONE(), vq4=dirac_slash(q4,4)-mq4*dirac_ONE();
	
	ex s4=m2q1+m2q2+m2q3+m2q4-s2-s3;
    scalar_products sp;
    sp.add(q2, q3, (s4-m2q2-m2q3)/2); 
    sp.add(q4, q3, (s2-m2q4-m2q3)/2); 
    sp.add(q2, q4, (s3-m2q2-m2q4)/2); 
  
    sp.add(q2, q2, m2q2);
    sp.add(q3, q3, m2q3);
    sp.add(q4, q4, m2q4);
    
    sp.add(h1,h1,-1);
    sp.add(h2,h2,-1);
    sp.add(h3,h3,-1);
    sp.add(h4,h4,-1);
    
    sp.add(h2,q2,0);
    sp.add(h3,q3,0);
    sp.add(h4,q4,0);
     
    multivector<ex,3> v(0,2,2,2);
	v[0][0][0]=dirac_gammaL(); v[0][0][1]=dirac_gammaR();
	v[0][1][0]=dirac_gammaR(); v[0][1][1]=dirac_gammaL();
	v[1][0][0]=dirac_gamma(jmu)*dirac_gammaL(); v[1][0][1]=dirac_gamma(jmu)*dirac_gammaL();
	v[1][1][0]=dirac_gamma(jmu)*dirac_gammaR(); v[1][1][1]=dirac_gamma(jmu)*dirac_gammaR();
	
	multivector<ex,7> traces(0,2,2,2,2,2,2,2);
	for(uint i=0;i<2;i++)
	  for(uint j=0;j<2;j++)
	    for(uint k=0;k<2;k++)
		  for(uint l=0;l<2;l++)
		  for(uint m=0;m<2;m++)
				for(uint n=0;n<2;n++){
				ex vik=v[i][k][0];
				ex vim=v[i][m][0].subs(mu==nu);
				ex vjl=v[j][l][1].subs(mu==alpha);
				ex vjn=v[j][n][1].subs(mu==beta);
	
				traces[i][j][k][l][m][n][0]=dirac_trace(vq2*vik*vq1*vjl)*dirac_trace(vq3*vim*vq4*vjn);
				traces[i][j][k][l][m][n][1]=-dirac_trace(vq2*vik*vq1*vjl*vq3*vim*vq4*vjn);
			}
	
	vector<ex> prop(2,0);
	prop[0]=1;
	prop[1]=lorentz_g(imu,inu);
	
	multivector<ex,2> prop2(0,2,2);
	for(uint i=0;i<2;i++)
	  for(uint j=0;j<2;j++){
				prop2[i][j]=prop[i]*prop[j].subs(lst(mu==alpha,nu==beta));
			}
	
	//ofstream f("M2.dat");					
	for(uint i=0;i<2;i++)
	  for(uint j=0;j<2;j++)
	    for(uint k=0;k<2;k++)
			for(uint l=0;l<2;l++)
				for(uint m=0;m<2;m++)
				for(uint n=0;n<2;n++)
						for(uint o=0;o<2;o++)
					{
							//cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<endl;
							M2[i][j][k][l][m][n][o]=(traces[i][j][k][l][m][n][o]*prop2[i][j]).simplify_indexed(sp);
							//cout<<M2[i][j][k][l][m][n][o]<<endl<<endl;
							//f<<M2[i][j][k][l][m][n][o]<<endl;
						}
}


ex get_integral(const multivector<ex,4>& a, const vector<ex>& mass, const vector<int>& op, double m1, double m2, double m3, double m4) const{
	
	ex q10=(s2+mq1*mq1-mq2*mq2)/(2*sqrt(s2)), lq1l=sqrt(q10*q10-mq1*mq1);
    ex q30=(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
    ex q20=(mq1*mq1+mq2*mq2-s2)/(2*mq1), lq2l=sqrt(q20*q20-mq2*mq2);
    
    ex total=0;
    for(uint i=0;i<a.size();i++) if(!mass[i].is_zero())
	  for(uint j=0;j<a.size();j++)
		for(uint k=0;k<2;k++)
		for(uint l=0;l<2;l++)
		for(uint m=0;m<2;m++)
		for(uint n=0;n<2;n++)
		for(uint r=0;r<2;r++)
		for(uint s=0;s<2;s++){ 
			ex coup=a[i][r][k][m]*a[j][s][l][n].conjugate();
			if(!coup.is_zero()){
					//cout<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
					ex integrand=M2[op[i]][op[j]][k][l][m][n][(r+s)%2];
					integrand=expand(integral(s3, mq1*mq1+mq3*mq3-2*q10*q30-2*lq1l*lq3l, mq1*mq1+mq3*mq3-2*q10*q30+2*lq1l*lq3l, integrand)\
													.eval_integ()/lq1l/sqrt(s2)*lq2l/mq1/mq1);
					double mm2=m2, mm3=m3;
					if(l) {mm2=m3; mm3=m2;}
					double result=ex_to<numeric>(integral(s2,std::pow(mm3+m4,2),std::pow(m1-mm2,2),integrand.subs(lst(mq1 == m1, mq2 == mm2, mq3 == mm3, mq4 == m4))).evalf()).to_double()/std::pow(M_PI,3)/512;
				    ex partial=result*coup/(pow(mass[i],2)*pow(mass[j],2));
				    //cout<<partial<<endl;
				    total=total+partial;
			}
			}
	
	return total;
}
ex get_integral_symb(const multivector<ex,4>& a, const vector<ex>& mass, const vector<int>& op, ex m1) const{
	
	ex q10=(s2+m1*m1)/(2*sqrt(s2)), lq1l=(m1*m1-s2)/(2*sqrt(s2));
    ex q30=sqrt(s2)/2, lq3l=q30;
    ex q20=(m1*m1-s2)/(2*m1), lq2l=q20;
    
    ex total=0;
    for(uint i=0;i<a.size();i++) if(!mass[i].is_zero())
	  for(uint j=0;j<a.size();j++)
		for(uint k=0;k<2;k++)
		for(uint l=0;l<2;l++)
		for(uint m=0;m<2;m++)
		for(uint n=0;n<2;n++)
		for(uint r=0;r<2;r++)
		for(uint s=0;s<2;s++){
			ex coup=a[i][r][k][m]*a[j][s][l][n].conjugate();
			if(!coup.is_zero()){
					//cout<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
					ex integrand=M2[op[i]][op[j]][k][l][m][n][(r+s)%2].subs(lst(mq1 == m1, mq2 == 0, mq3 == 0, mq4 == 0));
					integrand=expand(integral(s3, m1*m1-2*q10*q30-2*lq1l*lq3l, m1*m1-2*q10*q30+2*lq1l*lq3l, integrand)\
													.eval_integ()/lq1l/sqrt(s2)*lq2l/m1/m1);
					//integrand=integrand.subs(lst(mq1 == m1, mq2 == 0, mq3 == 0, mq4 == 0));
					//integrand=integrand.subs(pow(m1,2)/4-s2/2+pow(s2/m1,2)/4==pow((m1-s2/m1)/2,2));
					
										
					double mm2=0, mm3=0, m4=0;
					//if(l) {mm2=m3; mm3=m2;}
					ex result=integral(s2,std::pow(mm3+m4,2),pow(m1-mm2,2),integrand).eval_integ()/pow(Pi,3)/512;
					
				    ex partial=result*coup/(pow(mass[i],2)*pow(mass[j],2));
				    total=total+partial;
			}
			}
	return total;
}
ex get_integral_meson(const multivector<ex,4>& a, const vector<ex>& mass, const vector<int>& op, ex mm, ex m1, ex m2, ex m3, ex m4) const{
	
	ex q10=(s2+mq1*mq1-mq2*mq2)/(2*sqrt(s2)), lq1l=sqrt(q10*q10-mq1*mq1);
    ex q30=(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
    ex total=0;
    for(uint i=0;i<a.size();i++)
	  for(uint j=0;j<a.size();j++)
		for(uint k=0;k<2;k++)
		for(uint l=0;l<2;l++)
		for(uint m=0;m<2;m++)
		for(uint n=0;n<2;n++){
			ex coup=a[i][0][k][m]*a[j][0][l][n].conjugate();
		if(!coup.is_zero()){
				//cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<" "<<endl;
					ex integrand=M22[op[i]][op[j]][k][l][m][n];
					//cout<<collect_common_factors(expand(a[i][0][k][m]))<<endl;
					//cout<<collect_common_factors(expand(a[j][0][l][n].conjugate()))<<endl;
					integrand=expand(integral(s3, mq1*mq1+mq3*mq3-2*q10*q30-2*lq1l*lq3l, mq1*mq1+mq3*mq3-2*q10*q30+2*lq1l*lq3l, integrand)/lq1l/s2);
					ex result=integrand.subs(lst(sqrt(s2) == mm, s2==mm*mm, mq1 == m1, mq2 == m2, mq3 == m3, mq4 == m4))/Pi/128;
					ex mi=mass[i];
					if(mi.is_zero()) mi=mm;
					ex mj=mass[j];
					if(mj.is_zero()) mj=mm;
				    ex partial=result*coup/(pow(mi,2)*pow(mj,2));
				    //cout<<i<<"-"<<op[i]<<" "<<j<<"-"<<op[j]<<" "<<a[i]*a[j].conjugate()/(pow(mass[i],2)*pow(mass[j],2))<<endl<<endl;
					
				    total=total+partial;
			}
		}
	
	return total;
}

ex get_integral_meson2(const multivector<ex,4>& a, const vector<ex>& mass, const vector<int>& op, ex mm, ex m1, ex m2, ex m3, ex m4) const{
	
	ex q10=(s2+mq1*mq1-mq2*mq2)/(2*sqrt(s2)), lq1l=sqrt(q10*q10-mq1*mq1);
    ex q30=(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
    ex total=0;
    for(uint i=0;i<a.size();i++)
	  for(uint j=0;j<a.size();j++)
		for(uint k=0;k<2;k++)
		for(uint l=0;l<2;l++)
		for(uint m=0;m<2;m++)
		for(uint n=0;n<2;n++){
			ex coup=a[i][0][k][m]*a[j][0][l][n].conjugate();
		if(!coup.is_zero()){
				//cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<" "<<endl;
					ex integrand=M22[op[i]][op[j]][k][l][m][n];
					//cout<<collect_common_factors(expand(a[i][0][k][m]))<<endl;
					//cout<<collect_common_factors(expand(a[j][0][l][n].conjugate()))<<endl;
					integrand=expand(integral(s3, mq1*mq1+mq3*mq3-2*q10*q30-2*lq1l*lq3l, mq1*mq1+mq3*mq3-2*q10*q30+2*lq1l*lq3l, integrand)/lq1l/s2);
					ex result=integrand.subs(lst(sqrt(s2) == mm, s2==mm*mm, mq1 == m1, mq2 == m2, mq3 == m3, mq4 == m4))/Pi/128;
					ex mi=mass[i];
					if(mi.is_zero()) mi=mm;
					ex mj=mass[j];
					if(mj.is_zero()) mj=mm;
				    ex partial=result*coup/(pow(mi,2)*pow(mj,2));
				    //cout<<i<<"-"<<op[i]<<" "<<j<<"-"<<op[j]<<" "<<a[i]*a[j].conjugate()/(pow(mass[i],2)*pow(mass[j],2))<<endl<<endl;
					
				    total=total+partial;
			}
		}
	
	return total;
}
/*
ex get_integral_meson(const multivector<ex,2>& a, const vector<ex>& mass, const vector<int>& op, double meson_mass, double m3, double m4) const{
	
	ex q10=(meson_mass)/(2), lq1l=sqrt(q10*q10-mq1*mq1);
    ex q30=(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
    ex q20=(mq1*mq1+mq2*mq2-s2)/(2*mq1), lq2l=sqrt(q20*q20-mq2*mq2);
    
    
    ex total=0;
    for(uint i=0;i<a.size();i++)
	  for(uint j=0;j<a.size();j++)
		for(uint k=0;k<2;k++)
			for(uint l=0;l<2;l++) if(!(a[i][k]*a[j][l]).is_zero()){
					//cout<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
					ex integrand=M2[op[i]/2][op[j]/2][op[i]%2][op[j]%2][(k+l)%2];
					integrand=expand(integral(s3, mq1*mq1+mq3*mq3-2*q10*q30-2*lq1l*lq3l, mq1*mq1+mq3*mq3-2*q10*q30+2*lq1l*lq3l, integrand)\
													.eval_integ()/lq1l/sqrt(s2)*lq2l/mq1/mq1);
					double mm2=m2, mm3=m3;
					if(l) {mm2=m3; mm3=m2;}
					double result=ex_to<numeric>(integral(s2,pow(mm3+m4,2),pow(m1-mm2,2),integrand.subs(lst(mq1 == m1, mq2 == mm2, mq3 == mm3, mq4 == m4))).evalf()).to_double()/pow(M_PI,3)/512;
				    ex partial=result*a[i][k]*a[j][l].conjugate()/(pow(mass[i],2)*pow(mass[j],2));
				    //cout<<partial<<endl;
				    total=total+partial;
			}
	
	return total;
  }
*/  
multivector<ex,7> M2;
multivector<ex,6> M22;
realsymbol s2, s3;
realsymbol mq1, mq2, mq3, mq4;

};



#endif
