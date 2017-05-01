#ifndef MODEL_H
#define MODEL_H

#define _USE_MATH_DEFINES
#include <iostream>
#include <set>
#include <ginac/ginac.h>
#include <cmath>
#include <complex>
#include <memory>
#include <tr1/memory>
#include "TRandom3.h"
#include "TMath.h"

//g++ teste.cpp -o teste -lcln -lginac
using namespace std;
//using namespace tr1;
using namespace GiNaC;
///A class containing the value and uncertainty of an experimental measure
class measure{
public:
	measure(double v=0,double e=0):value(v),error(e){}
	measure operator*(measure m2) const{
		const measure & m1=*this;
		return measure(m1.value*m2.value,sqrt(std::pow(m1.value*m2.error,2)+std::pow(m2.value*m1.error,2)));
		}
	measure operator/(measure m2) const{
		const measure & m1=*this;
		return measure(m1.value/m2.value,sqrt(std::pow(m1.value/std::pow(m2.value,2)*m2.error,2)+std::pow(m1.error/m2.value,2)));
		}
	double value;	
	double error;
	};
///A base class representing an experimental measure
class observable{
public:
 observable(): copies(1) {}
 virtual ~observable(){}

/**\param hipothesis the theoretical hypothesis
 * \return the logarithm of the probability of measuring what was measured,
 * assuming that the hypothesis is true
 */
 virtual double loglikelihood(double hipothesis) const = 0;
 virtual double error(double hipothesis) const = 0;
 
 //virtual int type() const =0;
 mutable uint copies;
};

///An experimental measure which is an upper limit on a parameter with a given Confidence Level
class limitedobs: public observable{
public:
/**\param limit upper limit on the parameter
 * \param m minimum possible value of the parameter
 * \param p 1-Confidence Level 
 * */
 //limitedobs(double limit, double cl=0.9, double m=0): s(fabs(limit-m)/1.282), min(m),lim(limit) {
//	 if(cl==0.95) s*=1.282/1.645;
  limitedobs(double limit, double cl=0.9, double m=0): s(fabs(limit-m)/(1.282+sqrt(M_PI_2))), min(m),lim(limit) {
	 if(cl==0.95) s*=(1.282+sqrt(M_PI_2))/(1.645+sqrt(M_PI_2));
	 }
 ~limitedobs(){}
 double loglikelihood(double hipothesis) const {
 	 double diff=(hipothesis-min-sqrt(M_PI_2)*s)/s;
 	 if(diff<0) diff=0;
	 return diff*diff/2;
	}
  double error(double hipothesis) const {
	double diff=(hipothesis-min )/s;	
	return diff;
	}
 
 double s,min,lim;
};

///An experimental measure of a parameter which is a mean value and a standard deviation
class gaussobs: public observable{
public:

/**\param mean mean value of the measure
 * \param sigma standard deviation of the measure
 * */
 gaussobs(measure v): m(v.value), s(v.error) {}
 gaussobs(double mean, double sigma): m(mean), s(mean*sigma) {}
 ~gaussobs(){}
 double loglikelihood(double hipothesis) const {
	double diff=(m-hipothesis)/s;	
	return diff*diff/2;
	}
  double error(double hipothesis) const {
	double diff=(hipothesis-m)/s;	
	return diff;
	}
 const double m, s;
};


///the same as gaussobs but with a different initializer, such that the uncertainty sigma is absolute
class gauss2obs: public observable{
public:

/**\param mean mean value of the measure
 * \param sigma standard deviation of the measure
 * */
 gauss2obs(measure v): m(v.value), s(v.error) {}
 gauss2obs(double mean, double sigma): m(mean), s(sigma) {}
 ~gauss2obs(){}
 double loglikelihood(double hipothesis) const {
	double diff=(m-hipothesis)/s;	
	return diff*diff/2;
	}
 double error(double hipothesis) const {
	double diff=(hipothesis-m)/s;	
	return diff;
	}
 
 double expected()const {return m;}
 //int type() const {return 1;}
 const double m, s;
};

///A parameter which will be fitted in the simulation
class freeparameter{
public:
 /**\param mi minimum possible value for the parameter
  * \param ma maximum possible value for the parameter
  * \param r random number generator
  * */
 freeparameter(double mi, double ma, TRandom3 * r,double ss=1e-2): min(mi), max(ma), value(mi+(ma-mi)*r->Rndm()), step((ma-mi)*ss) {}
 
 ///changes randomly the ::value of the parameter, the standard deviation is ::step 
 /**\param r random number generator
  * */
 void next(TRandom3 * r,double f=1){
	 //double x=r->Gaus()*step;
	 value+=r->Gaus()*step*f;
	 }
  ///checks if the value of the parameter is between ::min and ::max
  bool isvalid() const {
	return min<=value && value<=max;
  }
  ///probability distribution, to be used by the Markov Chain Monte Carlo simulation
  /**\return \f$(\frac{x-::value}{::step})^2\f$
   * */
  double dist(double x) const {
	return std::pow((x-value)/step,2);
  }
  ///minimum possible value for the parameter
  double min;
  ///maximum possible value for the parameter
  double max;
  ///value of the parameter
  double value;
  ///standard deviation of the random changes of ::value in \ref next(TRandom3 *)
  double step;
 
};

///A parameter which will be fitted in the simulation
class discreteparameter{
public:
 /**\param mi minimum possible value for the parameter
  * \param ma maximum possible value for the parameter
  * \param r random number generator
  * */
  discreteparameter(int mi, int ma, TRandom3 * r): min(mi), max(ma), value(mi+r->Integer(ma-mi+1)) {}
  ///minimum possible value for the parameter
  double min;
  ///maximum possible value for the parameter
  double max;
  ///value of the parameter
  double value;
};

///vector of parameters
class parameters: public vector< freeparameter >{
public:

//parameters(){}
///changes randomly the value of the parameters 
void next(TRandom3 * r, double f=1){
	for(iterator i=begin();i!=end();i++){
		i->next(r,f);
		}
		
	//for(uint i=0;i<discrete.size();i++){
		//discrete[i].next(r);
		//}
	}
	
///checks if all the values are between their minimums and maximums
bool isvalid() const{
	for(const_iterator i=begin();i!=end();i++){
		if(!i->isvalid()) return 0;
		}
	return 1;
	}	
///checks if this and another vector of parameters are within 1sigma of distance 
double dist(const parameters& p) const{
	double total=0;
	for(uint i=0;i<size();i++){
		total+=at(i).dist(p[i].value);
		}
	return sqrt(total/size());
}

void setvalues(const parameters& p){
	
	for(uint i=0;i<size();i++){
		at(i).value=p[i].value;
		}
	}

double area() const{
	float a=1;
	for(const_iterator i=begin();i!=end();i++){
		a*=i->step;
		}
	return a;
	}

double gausslikelihood(const parameters & p2) const{
	double l=1;
	for(uint i=0;i<size();i++){
		l*=TMath::Gaus((p2[i].value-at(i).value)/at(i).step);
	}
	return l;
	}
	
	lst p;
	vector<double> values;
	//vector<discreteparameters> discrete;
};

///Base class to do the calculus of a constraint to the model
class calcu{
	public:
/**\param hipothesis the theoretical hypothesis
 * \return the logarithm of the probability of measuring what was measured,
 * assuming that the hypothesis is true
 */
 virtual double operator()(const parameters & p) const=0; 
 //virtual int type() const =0;
};

///class to do the calculus of a constraint based on a GiNaC compiled expression
class calcuba:public calcu{
	public:
	calcuba(observable * ob, const FUNCP_CUBA & e0): calcu(), o(ob), e(e0){}
	
	double operator()(const parameters & p) const{
		 double ret=1000;
	 int pass=1;
	 
	/* try{
		 ret=ex_to<numeric>(e.subs(p.p,subs_options::no_pattern).evalf()).to_double();
	 }
	 catch(GiNaC::pole_error e){
	  pass=0;
	  cout<<"Pole error"<<endl;
	 }
	 catch(...){
	  cout<<"Other exception"<<endl;
	  exit(1);
	 }
	 */
	 int n=p.values.size(), m=1;
	 e(&n,&(p.values[0]),&m,&ret);
	 if(pass) ret=o->loglikelihood(ret);
	 else ret=1000;
	 
	 return ret;
    }
    	
	shared_ptr<observable> o;
	FUNCP_CUBA e;
	};
	

///class to do the calculus of a constraint based on a GiNaC symbolic expression
class calcuex:public calcu{
	public:
	calcuex(observable * ob, const ex & e0): calcu(), o(ob), e(e0){}
	~calcuex(){}
	
	double operator()(const parameters & p) const{
		 double ret=1000;
	 int pass=1;
	 try{
		 ret=ex_to<numeric>(e.subs(p.p,subs_options::no_pattern).evalf()).to_double();
	 }
	 catch(GiNaC::pole_error e){
	  pass=0;
	  cout<<"Pole error"<<endl;
	 }
	 catch(exception e){
	  cout<<e.what()<<endl;
	  }
	 catch(...){
	  cout<<"Other exception"<<endl;
	  exit(1);
	 }
	 if(pass) ret=o->loglikelihood(ret);
	 else ret=1000;
	 
	 return ret;
    }
    
    double error(const parameters & p) const{
		 double ret=1000;
	 int pass=1;
	 try{
		 cout<<e<<endl;
		 cout<<e.subs(p.p)<<endl;
		 
		 ret=ex_to<numeric>(e.subs(p.p,subs_options::no_pattern).evalf()).to_double();
	 }
	 catch(GiNaC::pole_error er){
	  pass=0;
	  cout<<"Pole error"<<endl;
	 }
	 catch(exception er){
	  pass=0;
	  
	  cout<<er.what()<<endl;
	  cout<<e.subs(p.p,subs_options::no_pattern).evalf()<<endl;
	  }
	 catch(...){
	  cout<<"Other exception"<<endl;
	  exit(1);
	 }
	 if(pass) ret=o->error(ret);
	 else ret=1000;
	 
	 return ret;
    }
    	
shared_ptr<observable> o;
ex e;
};
	
///theoretical expression for an experimental measure
class prediction{
public:
 prediction(observable * ob, const FUNCP_CUBA & e0): calculate(new calcuba(ob,e0)) {}
 prediction(observable * ob, const ex & e0): calculate(new calcuex(ob,e0)) {}
 prediction(calcu * c): calculate(c) {}
 prediction(const prediction& p): calculate(p.calculate) {}
 ~prediction(){}
 
 double loglikelihood(const parameters & p) const { return (*calculate)(p);}
 ///theoretical expression for the experimental measure
 
 shared_ptr<calcu> calculate;
 //vector<ex> es;
};


///Abstract class for a model
class Model: public vector< prediction > {
public:
  
 Model(): r(new TRandom3(0)){}
 virtual ~Model(){delete r;};
 virtual parameters getlist(const parameters & p) const = 0;
 virtual parameters generateparameters(int max=0) const = 0;
 virtual int veto(const parameters & p, int max=0) const {return !p.isvalid();}
 
 ///calculates the probability of getting all the experimental measures if the model describes the reality
/**\param p vector with the values of the free parameters
 */

double likelihood(const parameters & p, bool check=1, int max=0) const{
	if(veto(p,max) && check) return 0;
	double total=loglike(p,0);
	if(total<-1000) return 0;
	return exp(total);
	}
 
double loglike(const parameters & p, bool check=1, int max=0) const{
	if(veto(p,max) && check) return -1000;
	parameters pp(getlist(p));
	
	double total=0; 
	int n=0;
	for(const_iterator i=begin();i!=end();i++) {
	try{
		n++;
		total+=i->loglikelihood(pp);
				}
    catch(exception e){
	  cout<<n<<e.what()<<endl;
	  exit(1);
	 }
	//catch(...){
		//cout<<"DD "<<n<<endl;
		//exit(1);
		//}
}

	return -total;
	}

  TRandom3 * r;
};

#endif
