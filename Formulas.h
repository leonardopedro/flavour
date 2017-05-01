#ifndef Formulas_H
#define Formulas_H

/*
#include "widthcalc.h"

#include "TH2F.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include <iostream>

#include "Math/Polynomial.h"
#include "Math/Interpolator.h"
#include <complex>
#include <cmath>
*/
using namespace std;

namespace BGLmodels{


enum FType{tLepton,tQuark};
enum FIsospin{iUp,iDown};
enum FFlavour{fElectron,fMuon,fTau,fAny};
enum FCharge{cParticle,cAntiParticle};
enum FHelicity{hLeft,hRight,hAny};
enum BSpin{sScalar, sVector,sAny};

/**
* @brief a fermion properties 
*/	
class Fermion{
public:

    Fermion(FType t, FIsospin i, FFlavour f=fAny, FCharge p=cParticle, FHelicity h=hAny): type(t), isospin(i), flavour(f), particle(p), helicity(h){}
    
    FType type;
    FIsospin isospin;
    FFlavour flavour;
    FCharge particle;
    FHelicity helicity;
};

/**
* @brief a meson properties
*/
class Meson{
public:

    Meson(const Fermion& qq1, const Fermion& qq2, ex m, ex d): q1(qq1), q2(qq2), mass(m), decay_factor(d){}
    
    Fermion q1, q2;
    ex mass;
    ex decay_factor;
};

constexpr double M_GF=1.166371e-5;
constexpr double M_MZ=91.1876;
constexpr double M_MW=80.398;
constexpr double M_cos2=std::pow(M_MW/M_MZ,2);
constexpr double M_Mu[3]={2.4e-3,1.29,172.9};
constexpr double M_Md[3]={5.3e-3,95e-3,4.2};
constexpr double M_Ml[3]={0.510998910e-3,105.6583715e-3,1776.82e-3};

typedef std::complex<double> CD;
typedef std::array<CD,3> Vector3c;
typedef std::array<std::array<CD,3>,3> Matrix3c;

/**
* @brief a class to represent the mixing matrices VCKM and VPMNS
*/
class Matrixx: public Matrix3c{
public:
Matrixx(): Matrix3c(){}
Matrixx(const Matrix3c& a): Matrix3c(a){}
///constructs a unitary Matrixx in the standard form
//Matrixx(std::initializer_list<Vector3c> l):

Matrixx(double c12, double c13, double c23,double s12, double s13, double s23, CD e13,CD e13t):
	Matrix3c({
		Vector3c({c12*c13,s12*c13,s13*e13t}),
		Vector3c({-s12*c23-c12*s23*s13*e13,c12*c23-s12*s23*s13*e13,s23*c13}),
		Vector3c({s12*s23-c12*c23*s13*e13,-c12*s23-s12*c23*s13*e13,c13*c23})
		}){} 

Matrixx(double t12, double t13, double t23, double d13): 
	Matrixx(std::cos(t12),std::cos(t13),std::cos(t23),std::sin(t12),std::sin(t13),std::sin(t23),std::exp(CD(0,d13)),std::exp(CD(0,-d13))){}
///computes the hermitian conjugate of the Matrixx 
const Matrixx conjugate() const{
	Matrixx res;
	for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++)
			res[i][j]=std::conj((*this)[j][i]);
	return res;
    }
};

const Matrixx Vnl=Matrixx(33.6*M_PI/180,9.11*M_PI/180,40.4*M_PI/180,M_PI/4).conjugate();
const Matrixx Vud(13.04*M_PI/180,0.201*M_PI/180,2.38*M_PI/180,1.2);


/**
* @brief definition of the couplings for the different BGL models
*/
class Mixes{
public:

  Mixes(ex tanb0,ex cp0,int genL=2,int genQ=2,int lup=0,int qup=0,int mssm=0):M(Matrix(),2,2),N(Matrix(),2,2),VN(Matrix(),2,2),N_(Matrix(),2,2),VN_(Matrix(),2,2),cp(cp0)
  {
	  tanb=tanb0;
	  M[tLepton][iDown]=Matrix(possymbol("m_e"),possymbol("m_\\mu"),possymbol("m_\\tau"));
	  M[tQuark][iUp]=Matrix(possymbol("m_u"),possymbol("m_c"),possymbol("m_t"));
	  M[tQuark][iDown]=Matrix(possymbol("m_d"),possymbol("m_s"),possymbol("m_b"));
	  const char * ln[3]={"1","2","3"};
	  const char * ll[3]={"e","\\mu","\\tau"};
	  const char * lu[3]={"u","c","t"};
	  const char * ld[3]={"d","s","b"};

	V.push_back(Matrix("U",ln,ll));
	V.push_back(Matrix("V",lu,ld));
	
	
	int up[2];
	up[0]=lup;
	up[1]=qup;
	
	vector< Matrix > delta;
	
	  vector<int> gL(3,0);
	  gL[genL]=1;
	  //Leptons
	  delta.push_back(Matrix(gL[0],gL[1],gL[2]));
	  //Quarks
	  vector<int> gQ(3,0);
	  gQ[genQ]=1;
	  delta.push_back(Matrix(gQ[0],gQ[1],gQ[2]));
	  
      for(uint i=0;i<2;i++){
		if(mssm){
		//Nu
	    N_[i][0]=0;
        //Nd
        N_[i][1]=0;
        //VNd
        VN_[i][1]=Matrix(tanb)*V[i];
        //Nu*V
        VN_[i][0]=Matrix(1/tanb)*V[i];
        }
        else if(up[i]){
		//Nu
		N_[i][0]=Matrix(tanb)+Matrix(-tanb-1/tanb)*V[i]*delta[i]*V[i].conjugate();
        //Nd
        N_[i][1]=Matrix(tanb)+Matrix(-tanb-1/tanb)*delta[i];
        //VNd
        VN_[i][1]=Matrix(tanb)*V[i]+Matrix(-tanb-1/tanb)*V[i]*delta[i];
        //Nu*V
        VN_[i][0]=Matrix(-1)*(Matrix(tanb)*V[i]+Matrix(-tanb-1/tanb)*V[i]*delta[i]);
		}else{
		//Nu
		N_[i][0]=Matrix(tanb)+Matrix(-tanb-1/tanb)*delta[i];
        //Nd
        N_[i][1]=Matrix(tanb)+Matrix(-tanb-1/tanb)*V[i].conjugate()*delta[i]*V[i];
        //VNd
        VN_[i][1]=Matrix(tanb)*V[i]+Matrix(-tanb-1/tanb)*delta[i]*V[i];
        //Nu*V
        VN_[i][0]=Matrix(-1)*(Matrix(tanb)*V[i]+Matrix(-tanb-1/tanb)*delta[i]*V[i]);
		}
		
		N[i][0]=N_[i][0]*M[i][0];
		N[i][1]=N_[i][1]*M[i][1];
		VN[i][1]=VN_[i][1]*M[i][1];
		VN[i][0]=M[i][0]*VN_[i][0];
      }
      appendtolst(replacements);  
      
  
  }
  ex mass(const Fermion& f) const{return M[f.type][f.isospin][f.flavour][f.flavour];}
  double massnum(const Fermion& f) const{return ex_to<numeric>(M[f.type][f.isospin][f.flavour][f.flavour].subs(replacements)).to_double();}
  
  void appendtolst(lst & reps) const{//,vector<ex>& var_errors
	  reps.append(M[0][1][0][0]==0.510998910e-3);
	  reps.append(M[0][1][1][1]==105.6583715e-3);
	  reps.append(M[0][1][2][2]==1776.82e-3);
	  
	  reps.append(M[1][0][0][0]==2.4e-3);
	  reps.append(M[1][0][1][1]==1.29);
	  reps.append(M[1][0][2][2]==172.9);
	  
	  reps.append(M[1][1][0][0]==5.3e-3);
	  reps.append(M[1][1][1][1]==95e-3);
	  reps.append(M[1][1][2][2]==4.2);
	  
	  vector< Matrix > Vn;
	Vn.push_back(Matrix(33.6*M_PI/180,9.11*M_PI/180,40.4*M_PI/180,M_PI/4).conjugate());
	Vn.push_back(Matrix(13.04*M_PI/180,0.201*M_PI/180,2.38*M_PI/180,1.2));
	for(uint i=0; i<2;i++)
		for(uint j=0; j<3;j++)
			for(uint k=0; k<3;k++)
				reps.append(V[i][j][k]==Vn[i][j][k]);
	  }
	
  lst reps;
  vector< Matrix > V;
  multivector<Matrix,2> M;
  multivector<Matrix,2> N;
  multivector<Matrix,2> VN;
 
  multivector<Matrix,2> N_;
  multivector<Matrix,2> VN_;
  
  lst replacements;
  ex cp;
  ex tanb;
};


/**
* @brief calculus of the constraints coming from the oblique parameters 
*/
class calcuOblique:public calcu{
	public:
	
	calcuOblique(): c1(0.741),c2(0.671), g1(c1*0.02,0.0397), g2(c2*0.02,0.1579) {}
	double operator()(const parameters & p) const{
	   double y=p[1].value;
	   double z=p[2].value;
	   double w=p[3].value;
	
	   double TT=(F(y*y,z*z)-F(w*w,z*z)+F(y*y,w*w))/(16*M_PI*M_MW*M_MW*(1-M_cos2));   
	   double Sparam=(std::pow(1-2*M_cos2,2)*G(y*y,y*y,M_MZ*M_MZ)+G(z*z,w*w,M_MZ*M_MZ)+2*log(z*w/y/y))/24/M_PI;
	   
	   double T1=c1*TT-c2*Sparam;
       double T2=c2*TT+c1*Sparam;
       
	return g1.loglikelihood(T1)+g2.loglikelihood(T2);
    }
	double F(double x,double y) const{
	 if(x==y) return 0;
	 return (x+y)/2-x*y*log(x/y)/(x-y);
	}
    double f(double t,double r) const{
	 if(r==0) return 0;
	 if(r<0) return 2*sqrt(-r)*atan(sqrt(-r)/t);
	 return sqrt(r)*log(fabs((t-sqrt(r))/(t+sqrt(r))));
	}

    double lnxy_xy(double x, double y) const{
	 if(x==y) return 1/y;
	 return log(x/y)/(x-y);
    }
    double G(double x,double y,double z) const{
	 double t=x+y-z;
	 double r=std::pow(z,2)-2*z*(x+y)+std::pow(x-y,2);
	 return -16.0/3+5*(x+y)/z-2*std::pow((x-y)/z,2)+r/std::pow(z,3)*f(t,r)+\
			3/z*lnxy_xy(x,y)*(std::pow(x,2)+std::pow(y,2)+(x-y)/std::pow(z,2)*(-std::pow(x,2)+std::pow(y,2)+std::pow(x-y,3)/3));
	}
   const double  c1,c2;
   const gauss2obs g1,g2;
};

constexpr double mt_mt=163.3,mt_mW=174.2,mt_mb=261.8;  

constexpr double C7SM(double x){
	return ((1/(x-1)+3)*x*std::log(x)+(-8*x*x-5*x+7)/6)*x/4/std::pow(x-1,3);
	}

constexpr double C8SM(double x){
	return (-3/(x-1)*x*std::log(x)+(-x*x+5*x+2)/2)*x/4/std::pow(x-1,3);
	}

constexpr double C7SM_MW=C7SM(std::pow(mt_mW/M_MW,2)),C7SM_Mt=C7SM(std::pow(mt_mt/M_MW,2)),C7SM_Mb=-0.353; 
constexpr double C8SM_MW=C8SM(std::pow(mt_mW/M_MW,2)),C8SM_Mt=C8SM(std::pow(mt_mt/M_MW,2)),C8SM_Mb=C8SM(std::pow(mt_mb/M_MW,2));


/**
* @brief calculus of the constraints coming from the b->s gamma decay 
*/
class calcubtosgamma2:public calcu{
	public:
	
constexpr static double calN=2.567e-3;  
constexpr static double a=7.8221,aee=0.4384,aer=-1.6981,a77=0.8161,a7r=4.8802,a7er=-0.7827,a88=0.0197,a8r=0.5680;
constexpr static double a8er=-0.0601,a87r=0.1923,a7i=0.3546,a8i=-0.0987,aei=2.4997,a87i=-0.0487,a7ei=-0.9067,a8ei=-0.0661;

	calcubtosgamma2(const Mixes& mixes): 
	  ii(2),
	  g1(3.43e-4,sqrt(2)*0.23e-4),
	  g2(9.2e-6,4e-6),
	  ratio(0){
		//cout<<"C7 "<<C7SM_Mt<<" "<<C7SM_MW<<" "<<C7SM(std::pow(261.8/M_MW,2))<<endl;
		double res[2];
	   constexpr double C7SM_[2]={C7SM_Mt,C7SM_Mt};
	   constexpr double C8SM_[2]={C8SM_Mt,C8SM_Mt};
	   for(uint j=0; j<2; j++){
		const uint i=2;
		const CD epsilon=conj(Vud[0][j])*Vud[0][i]/conj(Vud[2][j])/Vud[2][i];
		const double upsilon=norm(conj(Vud[2][j])*Vud[2][i]/Vud[1][i]);
		const CD R7=(C7SM_Mt)/C7SM_MW;
		const CD R8=(C8SM_Mt)/C8SM_MW;
		const CD R7_=0;
		const CD R8_=0;
		
		res[j]=a+aee*norm(epsilon)+aer*epsilon.real()+aei*epsilon.imag();
		res[j]+=a77*(norm(R7)+norm(R7_))+a7r*R7.real()+a7i*R7.imag();
		res[j]+=a88*(norm(R8)+norm(R8_))+a8r*R8.real()+a8i*R8.imag();
		res[j]+=a87r*(R8*conj(R7)+R8_*conj(R7_)).real()+a7er*(R7*conj(epsilon)).real()+a8er*(R8*conj(epsilon)).real();
		res[j]+=a87i*(R8*conj(R7)+R8_*conj(R7_)).imag()+a7ei*(R7*conj(epsilon)).imag()+a8er*(R8*conj(epsilon)).imag();
		res[j]*=calN/100*upsilon;
	}
	//cout<<"Btosgamma "<<res[0]/9.2e-6<<" "<<res[1]/3.15e-4<<endl;
		ifstream finter("interpolation.dat");
    
    if(!finter.is_open()){
		cout<<"ERROR: interpolation.dat not found"<<endl;
		exit(1);
		}
    vector<double> vinter0, vinter1, vinter2;
    while(!finter.eof()){
		double a=0,b=0,c=0;
		finter>>a>>b>>c;
		if(a!=0){
		// cout<<a<<" "<<b<<" "<<c<<endl;
		 vinter0.push_back(a);
		 vinter1.push_back(b);
		 vinter2.push_back(c);
	    }
		}
		
	inter1.SetData(vinter0,vinter1);
	inter2.SetData(vinter0,vinter2);
	
	finter.close();
	
	    ifstream finter2("masses.dat");
    
    if(!finter2.is_open()){
		cout<<"ERROR: masses.dat not found"<<endl;
		exit(1);
		}
    vector<vector<double> > m_(7);
    while(!finter2.eof()){
		for(uint i=0; i<7;i++) {
			double a=0;
			finter2>>a;
			
			if(a!=0) {
				if(i==0) a=log(a);
				else if(i<4) a*=1e-3;
				m_[i].push_back(a);
			//	cout<<a<<" ";
			}
		} //cout<<endl;
		}
	for(uint i=0; i<3;i++) {
		Md_[i].SetData(m_[0],m_[2*i+1]);
		Mu_[i].SetData(m_[0],m_[2*i+2]);
	}
//	cout<<"Eval "<<Mu_[2].Eval(log(100.0))<<endl;
//	cout<<"Eval "<<Md_[2].Eval(log(100.0))<<endl;
	
	finter2.close();

    ifstream finter3("interpolation2.dat");
    
    if(!finter3.is_open()){
		cout<<"ERROR: interpolation2.dat not found"<<endl;
		exit(1);
		}
    vector<double> vinter20, vinter21, vinter22;
    while(!finter3.eof()){
		double a=0,b=0,c=0;
		finter3>>a>>b>>c;
		if(a!=0){
		// cout<<a<<" "<<b<<" "<<c<<endl;
		 vinter20.push_back(a);
		 vinter21.push_back(b);
		 vinter22.push_back(c);
	    }
		}
		
	inter3.SetData(vinter20,vinter21);
	inter4.SetData(vinter20,vinter22);
	
	finter3.close();
	
	vector<ex> vex(24);
	
	   const uint i=ii; 
	for(uint j=0;j<2;j++)
	for(uint k=0;k<3;k++){
	    vex[j*6+k*2+0]=mixes.VN_[tQuark][iUp][k][j].conjugate()*mixes.VN_[tQuark][iUp][k][i];
	    vex[j*6+k*2+1]=mixes.N_[tQuark][iDown][j][k]*mixes.N_[tQuark][iDown][i][k].conjugate();
	    vex[j+k*2+12]=-mixes.VN_[tQuark][iUp][k][j].conjugate()*mixes.VN_[tQuark][iDown][k][i];
	    vex[j+k*2+18]=mixes.VN_[tQuark][iDown][k][j].conjugate()*mixes.VN_[tQuark][iDown][k][i];
		}
	lst l;	
    for(uint k=0;k<vex.size();k++){
		vex[k]=vex[k].subs(mixes.replacements).evalf();
		l.append(vex[k].real_part());
		l.append(vex[k].imag_part());
		}
	compile_ex(l, lst(mixes.tanb), fp);
	}
	
	double operator()(const parameters & p) const{
	 double tanb=p[0].value;
	 double y=p[1].value;
	   double z=p[2].value;
	   double w=p[3].value;
	double McH=y, MR=z, MI=w;
	
	double y0=y;
	if(y<mt_mt) y0=mt_mt;
	 double QCD1[2]={inter3.Eval(y0),inter1.Eval(y)};
	 double QCD2[2]={inter4.Eval(y0),inter2.Eval(y)};
	 
	 double Mu[3],Md[3];

	for(uint i=0;i<3;i++){
		Mu[i]=Mu_[i].Eval(log(y));
		Md[i]=Md_[i].Eval(log(z));
	}
	   const uint i=ii;
	   CD CC7[2],DD7[2],CC8[2],DD8[2];
	   double res[2];
	  // constexpr double C7SM_[2]={C7SM_Mt,C7SM_Mb};
	  // constexpr double C8SM_[2]={C8SM_Mt,C8SM_Mb};
	   
	   std::array<double,48> ret;
	   const int n=1,m=48;
	   fp(&n,&(tanb),&m,&(ret[0]));
	   for(uint j=0;j<2;j++){
	   	const double mbottom=Md[i];
		const double mstrange=Md[j];
			//ex mbottom=mixes.M[tQuark][iDown][i][i];
			//ex mstrange=mixes.M[tQuark][iDown][j][j];
	
			CD C7,D7,C8,D8;
			for(uint k=0;k<3;k++){
			double mup=Mu[k];
			double mdown=Md[k];
			//ex mup=mixes.M[tQuark][iUp][k][k];
			//ex mdown=mixes.M[tQuark][iDown][k][k];
			//f1+=
			double mmu=std::pow(mup/McH,2);
			double mmdR=std::pow(mdown/MR,2);
			double mmdI=std::pow(mdown/MI,2);
			
			double A0u=A0(mmu);
			double A1u=A1(mmu);
			double A2u=A2(mmu);
			double A3u=A3(mmu);
            double A0d=(A0(mmdR)+A0(mmdI));
			double A1d=(A1(mmdR)-A1(mmdI));
			
			CD f1(ret[j*12+4*k+0],ret[j*12+4*k+1]);
			C7+=f1*A2u;
			C8+=-2.0*f1*A0u;
			
			CD f2=CD(ret[36+j*2+4*k+0],ret[36+j*2+4*k+1])*mstrange*mbottom/mup/mup;
			//CD f2=f1*mstrange*mbottom/mup/mup;
			D7+=f2*A2u;
			D8+=-2.0*f2*A0u;
			
			CD f12(ret[24+j*2+4*k+0],ret[24+j*2+4*k+1]);			
			C7+=f12*A3u;
			C8+=2.0*f12*A1u;
			
			CD f4(ret[j*12+4*k+2],ret[j*12+4*k+3]);
			C7+=f4*A0d/3.0;
			C8+=-f4*A0d;
			
			C7+=f4*A1d/3.0;
			C8+=-f4*A1d;
			
			CD f6=f4*mstrange*mbottom/mdown/mdown;
			D7+=f6*A0d/3.0;
			D8+=-f6*A0d;
			}
		uint j0=j;	
		CC7[j]=(QCD1[j]*C7+QCD2[j]*C8)/2.0/conj(Vud[2][j])/Vud[2][i];
		DD7[j]=(QCD1[j]*D7+QCD2[j]*D8)/2.0/conj(Vud[2][j])/Vud[2][i];
		const double QCD3=(3*QCD2[j]/8+QCD1[j]);
		CC8[j]=QCD3*C8/2.0/conj(Vud[2][j])/Vud[2][i];
		DD8[j]=QCD3*D8/2.0/conj(Vud[2][j])/Vud[2][i];
		const CD epsilon=conj(Vud[0][j])*Vud[0][i]/conj(Vud[2][j])/Vud[2][i];
		const double upsilon=norm(conj(Vud[2][j])*Vud[2][i]/Vud[1][i]);
		const CD R7=(C7SM_Mt+CC7[j])/C7SM_MW;
		const CD R8=(C8SM_Mt+CD(0)*CC8[j])/C8SM_MW;
		const CD R7_=(DD7[j])/C7SM_MW;
		const CD R8_=CD(0)*(DD8[j])/C8SM_MW;
		
		
		res[j]=a+aee*norm(epsilon)+aer*epsilon.real()+aei*epsilon.imag();
		res[j]+=a77*(norm(R7)+norm(R7_))+a7r*R7.real()+a7i*R7.imag();
		res[j]+=a88*(norm(R8)+norm(R8_))+a8r*R8.real()+a8i*R8.imag();
		res[j]+=a87r*(R8*conj(R7)+R8_*conj(R7_)).real()+a7er*(R7*conj(epsilon)).real()+a8er*(R8*conj(epsilon)).real();
		res[j]+=a87i*(R8*conj(R7)+R8_*conj(R7_)).imag()+a7ei*(R7*conj(epsilon)).imag()+a8er*(R8*conj(epsilon)).imag();
		res[j]*=calN/100*upsilon;
		
		/*res[j]=a+aee*norm(epsilon)+aer*epsilon.real()+aei*epsilon.imag();
		res[j]+=a77*(norm(R7)+norm(R7_))+a7r*1+a7i*0;
		res[j]+=a88*(norm(R8)+norm(R8_))+a8r*R8.real()+a8i*R8.imag();
		res[j]+=a87r*1+a7er*(conj(epsilon)).real()+a8er*(R8*conj(epsilon)).real();
		res[j]+=a87i*0+a7ei*(conj(epsilon)).imag()+a8er*(R8*conj(epsilon)).imag();
		res[j]*=calN/100*upsilon;
		*/
	}		
	double r1=3.15e-4+0.00247*(norm(CC7[1])+norm(DD7[1])-0.706*CC7[1].real());
	
	//ratio=res[0]/9.2e-6;
	//cout<<"RATIO "<<ratio<<endl;
	return g1.loglikelihood(r1)+0*g2.loglikelihood(res[0]);
	   }
	
	double width(const parameters & p, int option=0) const{
	 double tanb=p[0].value;
	 double y=p[1].value;
	   double z=p[2].value;
	   double w=p[3].value;
	double McH=y, MR=z, MI=w;
	
	double y0=y;
	if(y<mt_mt) y0=mt_mt;
	 double QCD1[2]={inter3.Eval(y0),inter1.Eval(y)};
	 double QCD2[2]={inter4.Eval(y0),inter2.Eval(y)};
	 
	 double Mu[3],Md[3];

	for(uint i=0;i<3;i++){
		Mu[i]=Mu_[i].Eval(log(y));
		Md[i]=Md_[i].Eval(log(z));
	}
	   const uint i=ii;
	   CD CC7[2],DD7[2],CC8[2],DD8[2];
	   double res[2];
	  // constexpr double C7SM_[2]={C7SM_Mt,C7SM_Mb};
	  // constexpr double C8SM_[2]={C8SM_Mt,C8SM_Mb};
	   
	   std::array<double,24> ret;
	   const int n=1,m=24;
	   fp(&n,&(tanb),&m,&(ret[0]));
	   for(uint j=0;j<2;j++){
	   	const double mbottom=Md[i];
		const double mstrange=Md[j];
			//ex mbottom=mixes.M[tQuark][iDown][i][i];
			//ex mstrange=mixes.M[tQuark][iDown][j][j];
	
			CD C7,D7,C8,D8;
			for(uint k=0;k<3;k++){
			double mup=Mu[k];
			double mdown=Md[k];
			//ex mup=mixes.M[tQuark][iUp][k][k];
			//ex mdown=mixes.M[tQuark][iDown][k][k];
			//f1+=
			double mmu=std::pow(mup/McH,2);
			double mmdR=std::pow(mdown/MR,2);
			double mmdI=std::pow(mdown/MI,2);
			double A0u=0,A1u=0, A2u=0, A3u=0, A0d=0, A1d=0;
			
			if(option==0 || option==1){
			A0u=A0(mmu);
			A1u=A1(mmu);
			A2u=A2(mmu);
			A3u=A3(mmu);
			}
            if(option==0 || option==2){
            A0d=(A0(mmdR)+A0(mmdI));
			A1d=(A1(mmdR)-A1(mmdI));
			}
			if(option==3){
            A0d=(A0(mmdR));
			A1d=(A1(mmdR));
			}
			if(option==4){
            A0d=(A0(mmdI));
			A1d=(-A1(mmdI));
			}
			
			CD f1(ret[j*12+4*k+0],ret[j*12+4*k+1]);
			C7+=f1*A2u;
			C8+=-2.0*f1*A0u;
			
			CD f2=f1*mstrange*mbottom/mup/mup;
			D7+=f2*A2u;
			D8+=-2.0*f2*A0u;
			
			C7+=-f1*A3u;
			C8+=-2.0*f1*A1u;
			
			CD f4(ret[j*12+4*k+2],ret[j*12+4*k+3]);
			C7+=f4*A0d/3.0;
			C8+=-f4*A0d;
			
			C7+=f4*A1d/3.0;
			C8+=-f4*A1d;
			
			CD f6=f4*mstrange*mbottom/mdown/mdown;
			D7+=f6*A0d/3.0;
			D8+=-f6*A0d;
			
			}
		uint j0=j;	
		CC7[j]=(QCD1[j]*C7+QCD2[j]*C8)/2.0/conj(Vud[2][j])/Vud[2][i];
		DD7[j]=(QCD1[j]*D7+QCD2[j]*D8)/2.0/conj(Vud[2][j])/Vud[2][i];
		const double QCD3=(3*QCD2[j]/8+QCD1[j]);
		CC8[j]=QCD3*C8/2.0/conj(Vud[2][j])/Vud[2][i];
		DD8[j]=QCD3*D8/2.0/conj(Vud[2][j])/Vud[2][i];
		const CD epsilon=conj(Vud[0][j])*Vud[0][i]/conj(Vud[2][j])/Vud[2][i];
		const double upsilon=norm(conj(Vud[2][j])*Vud[2][i]/Vud[1][i]);
		const CD R7=(C7SM_Mt+CC7[j])/C7SM_MW;
		const CD R8=(C8SM_Mt+CD(0)*CC8[j])/C8SM_MW;
		const CD R7_=(DD7[j])/C7SM_MW;
		const CD R8_=CD(0)*(DD8[j])/C8SM_MW;
		
		
		res[j]=a+aee*norm(epsilon)+aer*epsilon.real()+aei*epsilon.imag();
		res[j]+=a77*(norm(R7)+norm(R7_))+a7r*R7.real()+a7i*R7.imag();
		res[j]+=a88*(norm(R8)+norm(R8_))+a8r*R8.real()+a8i*R8.imag();
		res[j]+=a87r*(R8*conj(R7)+R8_*conj(R7_)).real()+a7er*(R7*conj(epsilon)).real()+a8er*(R8*conj(epsilon)).real();
		res[j]+=a87i*(R8*conj(R7)+R8_*conj(R7_)).imag()+a7ei*(R7*conj(epsilon)).imag()+a8er*(R8*conj(epsilon)).imag();
		res[j]*=calN/100*upsilon;
		
		/*res[j]=a+aee*norm(epsilon)+aer*epsilon.real()+aei*epsilon.imag();
		res[j]+=a77*(norm(R7)+norm(R7_))+a7r*1+a7i*0;
		res[j]+=a88*(norm(R8)+norm(R8_))+a8r*R8.real()+a8i*R8.imag();
		res[j]+=a87r*1+a7er*(conj(epsilon)).real()+a8er*(R8*conj(epsilon)).real();
		res[j]+=a87i*0+a7ei*(conj(epsilon)).imag()+a8er*(R8*conj(epsilon)).imag();
		res[j]*=calN/100*upsilon;
		*/
	}		
	double r1=3.15e-4+0.00247*(norm(CC7[1])+norm(DD7[1])-0.706*CC7[1].real());
	
	//ratio=res[0]/9.2e-6;
	//cout<<"RATIO "<<ratio<<endl;
	return g1.error(r1);
	   }
	   
double A0(double x) const{
	return x*(2+3*x-6*x*x+ x*x*x+6*x*std::log(x))/(24*std::pow(1-x,4));
	}

double A1(double x) const{
	return x*(-3+4*x-x*x-2*std::log(x))/(4*std::pow(1-x,3));
	}

double A2(double x) const{
	return x/(6*std::pow(1-x,3))*((-7+5*x+8*x*x)/6.0+x*std::log(x)/(1-x)*(-2+3*x));
	}

double A3(double x) const{
	return (-3+8*x-5*x*x+(6*x-4)*std::log(x))*x/(6*std::pow(1-x,3));
	}
	
ROOT::Math::Interpolator inter1, inter2, inter3, inter4;  
ROOT::Math::Interpolator Mu_[3],Md_[3];

const uint ii;
FUNCP_CUBA fp;
 const gauss2obs g1,g2;
mutable double ratio;

};




/**
* @brief calculus of the constraints coming from the B->mu mu decay
*/
class calcuBmumu:public calcu{
	public:
	calcuBmumu(const Mixes& mix, const Meson & m,const Fermion& f3, const Fermion& f4, observable *ob, const char * name):
		meson(m),ff3(f3),ff4(f4),
		o(ob),
		gSr("gSr"),gSi("gSi"),gPr("gPr"),gPi("gPi"),gAr("gAr"),gAi("gAi"), mixes(mix) {
			const ex Nq=mixes.N[tQuark][m.q2.isospin][m.q1.flavour][m.q2.flavour];
			const ex Nq_=mixes.N[tQuark][m.q2.isospin][m.q1.flavour][m.q2.flavour].conjugate();
			const ex Nl=mixes.N[tLepton][f3.isospin][f3.flavour][f4.flavour];
			const ex Nl_=mixes.N[tLepton][f3.isospin][f4.flavour][f3.flavour].conjugate();
			possymbol MR("MR"),MI("MI"),McH("McH");		
			ex MR2=MR*MR,MI2=MI*MI,McH2=McH*McH;
			
			ex cLL=Nq_*Nl_*(1/MR2-1/MI2);
			ex cLR=Nq_*Nl*(1/MR2+1/MI2);
			ex cRL=Nq*Nl_*(1/MR2+1/MI2);
			ex cRR=Nq_*Nl*(1/MR2-1/MI2);
			
			ex ggS=-(2*M_GF/sqrt(2)*(-cRL-cRR+cLL+cLR)/4).subs(mixes.replacements).evalf();
			ex ggP=-(2*M_GF/sqrt(2)*(+cRL-cRR-cLL+cLR)/4).subs(mixes.replacements).evalf();
			CD ggA=0;
			if(m.q2.isospin==iDown && m.q2.flavour==2 && f3.flavour==1 && f4.flavour==1){
				ggA=-conj(Vud[2][m.q2.flavour])*Vud[2][m.q1.flavour]*Y(std::pow(M_Mu[2]/M_MW,2));
				ggA+=-conj(Vud[1][m.q2.flavour])*Vud[1][m.q1.flavour]*Y(std::pow(M_Mu[1]/M_MW,2));
				ggA*=M_GF*M_GF*M_MW*M_MW/M_PI/M_PI/2;
				
			}
			//ex gggA=0;
			//ex gggA=ggA.real()+I*ggA.imag();
						
			ex width=collect_common_factors(mesondwtest().subs(lst(gAr==ggA.real(),gAi==ggA.imag(),gSr==ggS.real_part(),gSi==ggS.imag_part(),gPr==ggP.real_part(),gPi==ggP.imag_part())).subs(mixes.replacements).evalf().real_part());
	
			compile_ex(lst(width), lst(mixes.tanb,McH,MR,MI), fp);
			
		}
	double operator()(const parameters & p) const{
	   //double factor=std::pow(M_GF*M_MW,4)/8/std::pow(M_PI,5)*std::sqrt(MM*MM-4*M_Ml[1]*M_Ml[1])*M_Ml[1]*M_Ml[1];
	   int n=4,m=1;
	   double ret=0;
	   fp(&n,&(p.values[0]),&m,&ret);
       
	   return o->loglikelihood(ret);
	}

double obsvalue(const parameters & p) const{
	   //double factor=std::pow(M_GF*M_MW,4)/8/std::pow(M_PI,5)*std::sqrt(MM*MM-4*M_Ml[1]*M_Ml[1])*M_Ml[1]*M_Ml[1];
	   int n=4,m=1;
	   double ret=0;
	   fp(&n,&(p.values[0]),&m,&ret);
       
	   return ret;
	}

	double Y(double x) const{
		return 1.0113*x/8/(1-x)*(4-x+3*x*log(x)/(1-x));
		}
	
ex mesondwtest() const{
	const Fermion& f1(meson.q2), f2(meson.q1);
	ex mesonmass=meson.mass;
	
	Fermion f3=ff3, f4=ff4;
	realsymbol q3("q3"), q4("q4");
	ex s2=pow(mesonmass,2);
		
		ex v1=0, v2=0;
		ex mq1=mixes.mass(f1),mq2=mixes.mass(f2);
		ex mq3=mixes.mass(f3),mq4=mixes.mass(f4);
		
		ex m2q1=mq1*mq1, m2q2=mq2*mq2, m2q3=mq3*mq3, m2q4=mq4*mq4;
		scalar_products sp;
		sp.add(q4, q3, (s2-m2q4-m2q3)/2);
		sp.add(q3, q3, m2q3);
		sp.add(q4, q4, m2q4);
		ex q10=(s2+mq1*mq1-mq2*mq2)/(2*sqrt(s2)), lq1l=sqrt(q10*q10-mq1*mq1);
		ex q30=(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
		
		ex a;
		a=-(gSr+I*gSi)*s2/(mq1+mq2);
		v1=v1+a*dirac_ONE();	
		v2=v2+a.conjugate()*dirac_ONE();
		a=-(gPr+I*gPi)*s2/(mq1+mq2);
		v1=v1+a*dirac_gamma5();
		v2=v2-a.conjugate()*dirac_gamma5();
		ex sl=(dirac_slash(q3,4)+dirac_slash(q4,4));
		a=(gAr+I*gAi);
		v1=v1+a*sl*dirac_gamma5();
		v2=v2+a.conjugate()*sl*dirac_gamma5();
		
	ex vq3=dirac_slash(q3,4)+mq3*dirac_ONE(), vq4=dirac_slash(q4,4)-mq4*dirac_ONE();
	ex dt=dirac_trace(vq3*v1*vq4*v2).simplify_indexed(sp);
	ex result=expand(dt*4*lq3l/s2/Pi/32);

	lst ltest;
	//ltest.append(conjugate(gL)==pow(abs(gL),2)/gL);
	//ltest.append(conjugate(gR)==pow(abs(gR),2)/gR);
//	ltest.append(conjugate(gS)==pow(abs(gS),2)/gS);
//	ltest.append(conjugate(gP)==pow(abs(gP),2)/gP);
//	ltest.append(conjugate(gA)==pow(abs(gA),2)/gA);
	
	return pow(meson.decay_factor,2)*collect_common_factors(result.subs(ltest));
	//return expand(ret.subs(lst(exp(-I*wild())==1/exp(I*wild()),sin(wild())==sqrt(1-pow(cos(wild()),2)))));
}
   const Meson meson;
   const Fermion& ff3, ff4;
   shared_ptr<observable> o;
   const realsymbol gSr,gSi, gPr,gPi,gAr,gAi;
   const Mixes mixes;
   FUNCP_CUBA fp;
   
};



/*

#include "defs.h"


//constexpr double Nll(int i, int j){return UL?(i==j)*(tanb*M_Ml[i]+(-tanb-1/tanb)*M_Ml[i]*(i==GL):((i==j)*tanb*M_Ml[i]+(-tanb-1/tanb)*conj(Vnl(GL,i))*Vnl(GL,j)*M_Ml[j]);}

#if UL==1
   #define Nl(i,j) ((i==j)*(tanb*M_Ml[i]+(-tanb-1/tanb)*M_Ml[i]*(i==GL)))
   #define VnlNl(i,j) (Vnl[i][j])*M_Ml[j]*(tanb+(-tanb-1/tanb)*(j==GL)))
#else
   #define Nl(i,j) ((i==j)*tanb*M_Ml[i]+(-tanb-1/tanb)*conj(Vnl[GL][i])*Vnl[GL][j]*M_Ml[j])
   #define VnlNl(i,j) (Vnl[i][j]*M_Ml[j]*(tanb+(-tanb-1/tanb)*(i==GL)))
#endif

#if UQ==1
   #define Nd(i,j) (M_Md[j]*(i==j)*(tanb+(-tanb-1/tanb)*(j==GQ)))
   #define Nu(i,j) (M_Mu[j]*((i==j)*tanb+(-tanb-1/tanb)*Vud[i][GQ]*conj(Vud[j][GQ])))
   #define VudNd(i,j) (Vud[i][j]*M_Md[j]*(tanb+(-tanb-1/tanb)*(j==GQ)))
   #define NuVud(i,j) (-M_Mu[i]*Vud[i][j]*(tanb+(-tanb-1/tanb)*(j==GQ)))
#else
   #define Nd(i,j) (M_Md[j]*((i==j)*tanb+(-tanb-1/tanb)*conj(Vud[GQ][i])*Vud[GQ][j]))
   #define Nu(i,j) (M_Mu[j]*(i==j)*(tanb+(-tanb-1/tanb)*M_Mu[j]*(j==GQ)))
   #define VudNd(i,j) (Vud[i][j]*M_Md[j]*(tanb+(-tanb-1/tanb)*(i==GQ)))
   #define NuVud(i,j) (-M_Mu[i]*Vud[i][j]*(tanb+(-tanb-1/tanb)*(i==GQ)))
#endif


*/
}
#endif
