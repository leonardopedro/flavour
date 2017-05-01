#include "MCMC.h"
#include "Juca.h"
#include "TF2.h"
#include "TProfile3D.h"
#include "THStack.h"
#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"



#include <cln/cln.h>
#include <cln/float.h>


int main(){
	Digits=5;
	cln::cl_inhibit_floating_point_underflow=1;
	
   //Int_t MyPalette[100];
   Double_t r[]    = {1, 0};
   Double_t g[]    = {1, 0};
   Double_t b[]    = {1, 0};
   Double_t stop[] = {0., 1.0};
   TColor::CreateGradientColorTable(2, stop, r, g,b, 100);
   
	//TH1F * pdf1=new TH1F("pdf1","pdf1",npoints,10,500);
	
	
	//TGraph * chi2=new TGraph(npoints);

	uint npoints=50;
	
	
		double llmax=-20,gmax=0;
		
		Juca* m=new Juca();
	Proposal prop4(m);
	prop4.findPeaks();
	
	cout<<"gp "<<m->gaussprob(prop4.floatPeak.pr)<<endl;
	lst l=m->getlist(prop4.floatPeak.pr);
	for(uint i=0; i< m->size();i++){
		double mean=m->at(i).calculate(l);
		cout<<i<<" "<<mean<<" "<<sqrt(2*m->at(i).o->loglikelihood(mean))<<endl;
	}
	for(uint i=0; i< prop4.floatPeak.pr.size();i++){
		cout<<i<<" "<<prop4.floatPeak.pr[i].value<<endl;
	}
	delete m;
	
	return 0;
}

