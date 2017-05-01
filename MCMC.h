#ifndef MCMC_H
#define MCMC_H

#include "model.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"

using namespace std;

///A class containing the parameters of a maximum of the likelihood function
class Peak {
public:
	Peak(const Model * m,int maxx=0): model(m), pr(m->generateparameters(maxx)), max(maxx){
		lmax=model->likelihood(pr);
		}
	
	void findPeak(){
	   llmax=model->loglike(pr,1,max);
	   area=1;
	  uint fixed=1e2;
	  uint f=fixed;
	  double d=1;
	  //cout<<"f "<<f<<"llmax "<<llmax<<endl;
		for(uint i=1e5;i;i--){
			parameters p1(pr);
			p1.next(model->r,d);
			
			double l1=llmax;
			if(!model->veto(p1,max)){l1=model->loglike(p1,1,max);
				}
			if(l1>llmax){pr=p1;  llmax=l1; f=fixed;}
			else {f--; if(!f) {d/=100; f=fixed; if(d<1e-2) break;}}
		    }
	   cout<<"d "<<d<<"llmax "<<llmax<<endl;
	   if(llmax<-1000) lmax=0;
	   else{lmax=exp(llmax);
		area=pr.area();
	    larea=area*lmax;
	    }
	   
	}
	/*
	void findPeak(){
	   llmax=model->loglike(pr);
	   area=1;
	  uint fixed=1e2;
	  uint f=fixed;
	  parameters p1(pr);
			
		for(uint j=1e4;j;j--){
			for(uint i=0;i<pr.size(); i++){
			   double s=pr[i].step/1e3;
			   pr[i].value+=s;
			   if(pr[i].value>pr[i].max){ s*=-1; pr[i].value+=2*s;}
			   double x=(model->loglike(pr)-llmax)*pr[i].step/1e2;
			   pr[i].value-=s;
			   if(fabs(x)>pr[i].step) x*=pr[i].step/fabs(x);
			   p1[i].value=pr[i].value+x;
			   if(p1[i].value>pr[i].max) p1[i].value=pr[i].max;
			   else if(p1[i].value<pr[i].min) p1[i].value=pr[i].min;
		   }
		   if(pr.dist(p1)<1e-6) {pr.setvalues(p1); f=0; break;}
		   cout<<"Loglike "<<llmax<<" "<<pr.dist(p1)<<endl;
		   pr.setvalues(p1);
		   llmax=model->loglike(pr);
	    }
			
	   if(f || llmax<-20) lmax=0;
	   else{
		   lmax=exp(llmax);
		area=pr.area();
	    larea=area*lmax;
	    }
	   
	}
	
	double RosenBrock(const double *xx )
	{	parameters p=model->generateparameters();
		for(uint i=0; i<p.size();i++) p[i].value=xx[i];
		
		return -model->loglike(pr);
	}
	
	void findPeak3(){
	   llmax=model->loglike(pr);
	   area=1;
	  uint fixed=1e2;
	  uint f=fixed;
	  parameters p1(pr);
	  ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS );
 
      min.SetMaxFunctionCalls(1000000);
      min.SetMaxIterations(100000);
      min.SetTolerance(0.001);
 
      ROOT::Math::Functor f(&RosenBrock,pr.size()); 
      double step[2] = {0.01,0.01};
      double variable[2] = { -1.,1.2};
	  char s[3]="x0";
	  
      min.SetFunction(f);
	fo
      // Set the free variables to be minimized!
      min.SetVariable(0,"x",variable[0], step[0]);
      min.SetVariable(1,"y",variable[1], step[1]);
 
   min.Minimize();
   	
		for(uint j=1e4;j;j--){
			for(uint i=0;i<pr.size(); i++){
			   double s=pr[i].step/1e3;
			   pr[i].value+=s;
			   if(pr[i].value>pr[i].max){ s*=-1; pr[i].value+=2*s;}
			   double x=(model->loglike(pr)-llmax)*pr[i].step/1e2;
			   pr[i].value-=s;
			   if(fabs(x)>pr[i].step) x*=pr[i].step/fabs(x);
			   p1[i].value=pr[i].value+x;
			   if(p1[i].value>pr[i].max) p1[i].value=pr[i].max;
			   else if(p1[i].value<pr[i].min) p1[i].value=pr[i].min;
		   }
		   if(pr.dist(p1)<1e-6) {pr.setvalues(p1); f=0; break;}
		   cout<<"Loglike "<<llmax<<" "<<pr.dist(p1)<<endl;
		   pr.setvalues(p1);
		   llmax=model->loglike(pr);
	    }
			
	   if(f || llmax<-20) lmax=0;
	   else{
		   lmax=exp(llmax);
		area=pr.area();
	    larea=area*lmax;
	    }
	   
	}
	*/
	bool adjuststeps(){
		 parameters p1(pr);
		 for(uint i=0;i<pr.size(); i++){
			   double s=p1[i].step;
			   p1[i].value+=s;
			   double x=(lmax-model->likelihood(p1))*2/lmax/s/s;
			   double x0=std::pow(2/(pr[i].max-pr[i].min),2);
			   if(x<x0) return 0;
			  // cout<<"X "<<x<<endl;
			   pr[i].step=1/sqrt(x);
			   p1[i].value-=s;
			   }
	   area=pr.area();
	   larea=area*lmax;
	   return 1;
	}
	
	const Model * model;
	parameters pr;
	double lmax, llmax;
	double area;
	double larea;
	bool max;
};

///A class containing the parameters of a proposal for the next step in the Markov Chain
class Proposal{
public:

Proposal(const Model * m): model(m), floatPeak(m),proposal(m){}

void findPeaks(uint ns=1, int max=0) {
	//float pmin=-100, pmax=100, s=0.1;
    //floatPeak.s=s;
    //floatPeak.lmax=0;
    //int imax=-1;
    floatPeak=Peak(model,max);
    floatPeak.lmax=0;
    floatPeak.llmax=-1000;
   cout<<"started"<<endl;
	//for(uint i=5e1;i;i--){
	for(uint i=ns;i;i--){
		Peak pp(model,max);
		pp.findPeak();
		if(pp.llmax>-15){
		//for(uint j=0; j< pp.pr.size();j++){
		//cout<<j<<" "<<pp.pr[j].value<<endl;
		//}
		//lst l=model->getlist(pp.pr);
		//for(uint j=0; j< model->size();j++){
		//	double mean=model->at(j).calculate(l);
		//cout<<j<<" "<<mean<<" "<<sqrt(2*model->at(j).o->loglikelihood(mean))<<endl;
		//}
		}
		if(pp.lmax>floatPeak.lmax){
			cout<<i<<" "<<pp.lmax<<endl;
			floatPeak.lmax=pp.lmax;
			floatPeak.pr=pp.pr;
			}
	}
	floatPeak.area=floatPeak.pr.area();
	
}

void getProposal(){
	if(model->r->Rndm()<=0.9) { 
		proposal.pr=floatPeak.pr;
		proposal.pr.next(model->r);
		return;
		}
		
	proposal.pr=model->generateparameters();
}

void getNextPoint(){
	getProposal();
	double l1=0;
	l1=model->likelihood(proposal.pr);
	if(model->r->Rndm()<=l1/floatPeak.lmax){
		floatPeak.lmax=l1;
		floatPeak.pr.setvalues(proposal.pr);
	}
}

const Model * model;

vector<Peak> vPeak;
Peak floatPeak, proposal;
double total;
};

#endif
