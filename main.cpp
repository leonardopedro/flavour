/**
 * @file main.c
 * @author Leonardo Pedro
 * @date May 2014
 * @brief Main file
 *
 * Here typically goes a more extensive explanation of what the header
 * defines. Doxygens tags are words preceeded by either a backslash @\
 * or by an at symbol @@.
 * @see http://www.stack.nl/~dimitri/doxygen/docblocks.html
 * @see http://www.stack.nl/~dimitri/doxygen/commands.html
 */

#include "MCMC.h"
#include "BGL.h"
#include "TF2.h"
#include "TProfile3D.h"
#include "THStack.h"
#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TFile.h"
#include "TVector.h"

#include <cln/cln.h>
#include <cln/float.h>

using namespace BGLmodels;

/** \mainpage Introduction
* The program produces figures presenting 68%, 95% and 99% CL allowed regions in parameter space. 
* To wit, we represent regions where the specific BGL model is able to fit the imposed experimental information at least as well as the corresponding goodness levels. 
* Some comments are in order. 
* This procedure corresponds to the profile likelihood method. 
* In brief, for a model with parameters \f$\vec p\f$, we compute the predictions for the considered set of observables \f$\vec O_{\mathrm{Th}}(\vec p)\f$. 
* Then, using the experimental information \f$\vec O_{\mathrm{Exp}}\f$ available for those observables, we build a likelihood function \f$\mathcal L(\vec O_{\mathrm{Exp}}|\vec O_{\mathrm{Th}}(\vec p))\f$ 
* which gives the probability of obtaining the experimental results \f$\vec O_{\mathrm{\mathrm{Exp}}}\f$ assuming that the model is correct.
* The likelihood function \f$\mathcal L(\vec O_{\mathrm{Exp}}|\vec O_{\mathrm{Th}}(\vec p))\f$ 
* encodes all the information on how the model is able to reproduce the observed data all over parameter space.
* Nevertheless, the knowledge of \f$\mathcal L(\vec O_{\mathrm{Exp}}|\vec O_{\mathrm{Th}}(\vec p))\f$ in a multidimensional parameter space
* can be hardly represented and one is led to the problem of reducing that information to one or two-dimensional subspaces. 
* In the profile likelihood method, for each point in the chosen subspace, the highest likelihood over the complementary, marginalized space, is retained. Let us clarify that likelihood 
* -- or chi-squared \f$\chi^2\equiv -2\log \mathcal L\f$ -- 
* profiles and derived regions such as the ones we represent, are thus insensitive to the size of the space over which one marginalizes; 
* this would not be the case in a Bayesian analysis, where an integration over the marginalized space is involved. The profile likelihood method seems adequate to our purpose,
* which is none other than exploring where in parameter space are the different BGL models able to satisfy experimental constraints,
* without weighting in eventual fine tunings of the models or parameter space volumes.
* For the numerical computations the libraries GiNaC and ROOT are used. *
*/


/**
* @brief the main function takes the arguments gL gQ lup qup which specify a BGL model and runs the simulation for that model, generating a ROOT file as the output 
*/
int main(int argc, char* argv[]){
    // Check the number of parameters
    if(argc < 5){
        std::cerr<<"Usage: "<<argv[0]<<" gL gQ lup qup [mssm]"<<std::endl;
        return 1;}
       
   
    int gL=atoi(argv[1]);
    int gQ=atoi(argv[2]);
    int lup=atoi(argv[3]);
    int qup=atoi(argv[4]);
    int mssm=0;
    if(argc>5) mssm=atoi(argv[5]);
    

   Digits=5;
   cln::cl_inhibit_floating_point_underflow=1;
   
   string ll[2][3]={{"#nu_{1}","#nu_{2}","#nu_{3}"},{"e","#mu","#tau"}};
   string qq[2][3]={{"u","c","t"},{"d","s","b"}};    
   //Int_t MyPalette[100];
   Double_t r[]    = {1, 0.3};
   Double_t g[]    = {1, 0.3};
   Double_t b[]    = {1, 0.3};
   Double_t stop[] = {0., 1.0};
   TColor::CreateGradientColorTable(2, stop, r, g,b, 100);
   //TH1F * pdf1=new TH1F("pdf1","pdf1",npoints,10,500);
   //TGraph * chi2=new TGraph(npoints);
   
    uint npoints=200;
    double init1=-3, final1=3;
    double init2=10, final2=1000;
    double initBmumu=0, finalBmumu=2;
    double initBsmumu=0, finalBsmumu=5;
    
	multivector<BGL *,4> ms(0,3,3,2,2);
	multivector<parameters,2> plots(parameters(),npoints,npoints);
	multivector<double,2> likely(-1000,npoints,npoints);
	multivector<double,2> likelyB(-1000,npoints,npoints);
	multivector<double,2> likely2(-1000,npoints,npoints);
	multivector<double,2> likely3(-1000,npoints,npoints);
	multivector<double,2> likely4(-1000,npoints,npoints);
	multivector<double,2> likely5(-1000,npoints,npoints);
	multivector<double,2> likely6(-1000,npoints,npoints);
	
	//ofstream mass("mass.out");
	
	//for(int gL=2;gL>=0;gL--)
	//for(int gQ=2;gQ>=0;gQ--)
	//for(int lup=0;lup<2;lup++)
	//for(int qup=0;qup<2;qup++){
		
//	for(uint qup=0;qup<2;qup++)
//for(uint gL=0;gL<3;gL++)
//for(uint lup=0;lup<2;lup++)
//for(uint gQ=0;gQ<3;gQ++){
		//if(gL==0 && gQ==2 && lup==0 && qup==1) {t=1; continue;}
		//if(t==0) continue;
	double llmax=-1000,gmax=0, McHmax=1000,MRmax=1000,MImax=1000,tbmax=1;
		BGL* m=new BGL(gL,gQ,lup,qup);
		m->mmmax=300;
		ms[gL][gQ][lup][qup]=m;
		
		char name[5]="0000";
		name[0]+=gL;
		name[1]+=gQ;
		name[2]+=lup;
		name[3]+=qup;
			
	/*
	TF1 * f1 = new TF1("f1",m,&BGL::BranchingRatio,-3,3,0,"BGL","BranchingRatio"); 
	
	TCanvas * c1=new TCanvas("c1","",800,600);
	f1->Draw();
	c1->SaveAs("BR.png");	
	TF2 * f2 = new TF2("f2",m,&BGL::topBranchingRatio,-3,3,80,175,0,"BGL","topBranchingRatio"); 
	
	TCanvas * c3=new TCanvas("c3","",800,600);
	f2->Draw("colz");
	c3->SaveAs("topBR.png");
	*/
	
	Proposal prop4(m);
	//m->stepsize=1e-4;
	prop4.findPeaks(100,1);
	llmax=log(prop4.floatPeak.lmax);
	 tbmax=prop4.floatPeak.pr[0].value;
	  McHmax=prop4.floatPeak.pr[1].value;
	  MRmax=prop4.floatPeak.pr[2].value;
	  MImax=prop4.floatPeak.pr[3].value;

	m->stepsize=1e-4;
	prop4.findPeaks(100);
	double llmax_=log(prop4.floatPeak.lmax);
	if(llmax_>llmax) {llmax=llmax_; 
	  tbmax=prop4.floatPeak.pr[0].value;
	  McHmax=prop4.floatPeak.pr[1].value;
	  MRmax=prop4.floatPeak.pr[2].value;
	  MImax=prop4.floatPeak.pr[3].value;
	}
	
	TH2F * limits4=new TH2F("limits4","Likelihood",npoints,init1,final1,npoints,init2,final2);
	TH2F * Bmumu_Bsmumu=new TH2F("Bmumu_Bsmumu","Likelihood",npoints,initBmumu,finalBmumu,npoints,initBsmumu,finalBsmumu);
	TH2F * limits_tb_MR=new TH2F("limits_tb_MR","Likelihood",npoints,init1,final1,npoints,init2,final2);
	TH2F * limits_tb_MI=new TH2F("limits_tb_MI","Likelihood",npoints,init1,final1,npoints,init2,final2);
	TH2F * limits_MR_MI=new TH2F("limits_MR_MI","Likelihood",npoints,init2,final2,npoints,init2,final2);
	TH2F * limits_MR_McH=new TH2F("limits_MR_McH","Likelihood",npoints,init2,final2,npoints,init2,final2);
	TH2F * limits_MI_McH=new TH2F("limits_MI_McH","Likelihood",npoints,init2,final2,npoints,init2,final2);
	
	
	for(uint i=0;i<npoints;i++)
		for(uint j=0;j<npoints;j++) {
			limits4->SetBinContent(i+1,j+1,-1000);
		    Bmumu_Bsmumu->SetBinContent(i+1,j+1,-1000);
		    limits_MR_McH->SetBinContent(i+1,j+1,-1000);
		    limits_MI_McH->SetBinContent(i+1,j+1,-1000);
			limits_MR_MI->SetBinContent(i+1,j+1,-1000);
			limits_tb_MR->SetBinContent(i+1,j+1,-1000);
			plots[i][j]=m->generateparameters();
			}
	uint steps=40e6;
	//uint steps=npoints*npoints;
	double brtaunu=1;
	for(uint i=steps;i;i--){
		//prop1.getNextPoint();
		//prop2.getNextPoint();
		//prop3.getNextPoint();
		double total=0;
		double gp=0;
		//if(i==steps) cout<<" total "<<m->loglike(prop4.floatPeak.pr)<<endl;
		//else{
		if(i>npoints*npoints){
			if(i==steps/2) {
				m->mmmax=1000;
				m->stepsize=1e-2;
				prop4.findPeaks(100);
				llmax_=log(prop4.floatPeak.lmax);
			         if(llmax_>llmax) {llmax=llmax_; 
				   tbmax=prop4.floatPeak.pr[0].value;
	  McHmax=prop4.floatPeak.pr[1].value;
	  MRmax=prop4.floatPeak.pr[2].value;
	  MImax=prop4.floatPeak.pr[3].value;
	}
			}
			if(i<steps) prop4.getNextPoint();
			total=log(prop4.floatPeak.lmax);
			//total=m->loglike(prop4.floatPeak.pr);
			//gp=m->gaussprob(prop4.floatPeak.pr);	
			}
		else{
			uint x=(i-1)/npoints, y=(i-1)%npoints;
			prop4.floatPeak.pr[0].value=init1+((x+0.5)*(final1-init1))/npoints;
			prop4.floatPeak.pr[1].value=init2+((y+0.5)*(final2-init2))/npoints;
			prop4.floatPeak.pr[2].value=0;
			prop4.floatPeak.pr[3].value=0;
			total=m->loglike(prop4.floatPeak.pr);
			}
		//}
	
		
		//gp=m->gaussprob_noT(prop4.floatPeak.pr);	
		//double Btaunu=m->loglike(prop1.floatPeak.pr,m->iBtaunu);
		//double BDtaunu=m->loglike(prop2.floatPeak.pr,m->iBDtaunu);
		//double BD2taunu=m->loglike(prop3.floatPeak.pr,m->iBD2taunu);
		double brtaunu0=ex_to<numeric>(m->BR_Htotaunu.subs(m->getlist(prop4.floatPeak.pr).p).evalf()).to_double();
		if(brtaunu0<brtaunu) brtaunu=brtaunu0;
		if(total>llmax) {llmax=total; gmax=gp;
		  tbmax=prop4.floatPeak.pr[0].value;
		  McHmax=prop4.floatPeak.pr[1].value;
		  MRmax=prop4.floatPeak.pr[2].value;
		  MImax=prop4.floatPeak.pr[3].value;
		}
		
		
		
		uint p1=uint((prop4.floatPeak.pr[0].value-init1)/(final1-init1)*npoints);
		double mr=prop4.floatPeak.pr[1].value;
		uint p2=uint((mr-init2)/(final2-init2)*npoints);
		mr+=prop4.floatPeak.pr[2].value;
		uint p3=uint((prop4.floatPeak.pr[2].value+prop4.floatPeak.pr[1].value-init2)/(final2-init2)*npoints);
		mr+=prop4.floatPeak.pr[3].value;
		uint p4=uint((mr-init2)/(final2-init2)*npoints);
		if(p1<npoints && p2<npoints && p3<npoints &&  p4<npoints){
		if(total>likely[p1][p2]){
			likely[p1][p2]=total;
			plots[p1][p2].setvalues(prop4.floatPeak.pr);
			limits4->SetBinContent(p1+1,p2+1,total);
		}
		if(total>likely2[p3][p4]){
			likely2[p3][p4]=total;
			//plots[gL][gQ][lup][qup][p1][p2].setvalues(prop4.floatPeak.pr);
			limits_MR_MI->SetBinContent(p3+1,p4+1,total);
		}
		
		if(total>likely3[p3][p2]){
			likely3[p3][p2]=total;
			//plots[gL][gQ][lup][qup][p1][p2].setvalues(prop4.floatPeak.pr);
			limits_MR_McH->SetBinContent(p3+1,p2+1,total);
		}
		
		if(total>likely4[p1][p3]){
			likely4[p1][p3]=total;
			//plots[gL][gQ][lup][qup][p1][p2].setvalues(prop4.floatPeak.pr);
			limits_tb_MR->SetBinContent(p1+1,p3+1,total);
		}
		if(total>likely5[p4][p2]){
			likely5[p4][p2]=total;
			//plots[gL][gQ][lup][qup][p1][p2].setvalues(prop4.floatPeak.pr);
			limits_MI_McH->SetBinContent(p4+1,p2+1,total);
		}
		parameters ppr=m->getlist(prop4.floatPeak.pr);
		
		double Bmumu=m->cBmumu->obsvalue(ppr)/(1e-10*m->planck/1.519e-12);
		double Bsmumu=m->cBsmumu->obsvalue(ppr)/(1e-9*m->planck/1.516e-12);
		//cout<<Bmumu<<" "<<Bsmumu<<endl;
		uint pB=uint((Bmumu-initBmumu)/(finalBmumu-initBmumu)*npoints);
		uint pBs=uint((Bsmumu-initBsmumu)/(finalBsmumu-initBsmumu)*npoints);
		if(pB<npoints && pBs<npoints)
		if(total>likelyB[pB][pBs]){
			likely[pB][pBs]=total;
			Bmumu_Bsmumu->SetBinContent(pB+1,pBs+1,total);
		}
		
	}else cout<<"OUT!!!"<<endl;
		
		if(i%100000==0){ 
			cout<<"Steps "<<i<<" logtb "<<prop4.floatPeak.pr[0].value<<" logMcH "<<prop4.floatPeak.pr[1].value;
			cout<<" total "<<total<<" gmax "<<gp<<" "<<brtaunu<<endl;
			}
		
		}
	//}
	
	TFile f((string("h")+string(name)+string(".root")).c_str(),"recreate"); 
    TVectorD v(5);
	v[0] = llmax; 
        v[1]=tbmax;
	v[2]=McHmax;
	v[3]=MRmax;
	v[4]=MImax;
    v.Write("vllmax");
    limits4->Write();
    Bmumu_Bsmumu->Write();
    limits_MR_MI->Write();
    limits_MR_McH->Write();
    limits_MI_McH->Write();
    limits_tb_MR->Write();
	f.Close();
	
	//for(int gL=2;gL>=0;gL--)
	//for(int gQ=2;gQ>=0;gQ--)
	//for(uint lup=0;lup<2;lup++)
	//for(uint qup=0;qup<2;qup++)
	uint min1=npoints, min2=npoints, min3=npoints;
	uint min11=npoints, min21=npoints, min31=npoints;
	uint min12=npoints, min22=npoints, min32=npoints;
	
	for(uint i=0;i<npoints;i++)
	for(uint j=0;j<npoints;j++){
			int binmax=limits4->GetBin(i+1,j+1);
			double rest=limits4->GetBinContent(binmax);
			if(rest>=llmax) rest=1;
			else rest=TMath::Prob(-2*(rest-llmax),2);
			if(rest>=0.05 && j<min1){min1=j;}
			if(rest>=0.05 && j<min11 && j>uint(180*npoints/990)){min11=j;}
			if(rest>=0.05 && j<min12 && j>uint(400*npoints/990)){min12=j;}
			limits4->SetBinContent(i+1,j+1,rest);
			
			rest=Bmumu_Bsmumu->GetBinContent(binmax);
			if(rest>=llmax) rest=1;
			else rest=TMath::Prob(-2*(rest-llmax),2);
			Bmumu_Bsmumu->SetBinContent(i+1,j+1,rest);			
			
			rest=limits_MR_MI->GetBinContent(binmax);
			if(rest>=llmax) rest=1;
			else rest=TMath::Prob(-2*(rest-llmax),2);
			limits_MR_MI->SetBinContent(i+1,j+1,rest);
			
			rest=limits_MR_McH->GetBinContent(binmax);
			if(rest>=llmax) rest=1;
			else rest=TMath::Prob(-2*(rest-llmax),2);
			if(rest>=0.05 && i<min2){min2=i;}
			if(rest>=0.05 && i<min21 && j>uint(180*npoints/990)){min21=i;}
			if(rest>=0.05 && i<min22 && j>uint(400*npoints/990)){min22=i;}
			limits_MR_McH->SetBinContent(i+1,j+1,rest);
			
			rest=limits_MI_McH->GetBinContent(binmax);
			if(rest>=llmax) rest=1;
			else rest=TMath::Prob(-2*(rest-llmax),2);
			if(rest>=0.05 && i<min3){min3=i;}
			if(rest>=0.05 && i<min31 && j>uint(180*npoints/990)){min31=i;}
			if(rest>=0.05 && i<min32 && j>uint(400*npoints/990)){min32=i;}
			limits_MI_McH->SetBinContent(i+1,j+1,rest);
			
			
			rest=limits_tb_MR->GetBinContent(binmax);
			if(rest>=llmax) rest=1;
			else rest=TMath::Prob(-2*(rest-llmax),2);
			limits_tb_MR->SetBinContent(i+1,j+1,rest);
		}
		double mmin1=init2+((min1+0.5)*(final2-init2))/npoints;
		double mmin2=init2+((min2+0.5)*(final2-init2))/npoints;
		double mmin3=init2+((min3+0.5)*(final2-init2))/npoints;
		
		double mmin11=init2+((min11+0.5)*(final2-init2))/npoints;
		double mmin21=init2+((min21+0.5)*(final2-init2))/npoints;
		double mmin31=init2+((min31+0.5)*(final2-init2))/npoints;
		
		double mmin12=init2+((min12+0.5)*(final2-init2))/npoints;
		double mmin22=init2+((min22+0.5)*(final2-init2))/npoints;
		double mmin32=init2+((min32+0.5)*(final2-init2))/npoints;
		
		//mass<<name<<" "<<mmin1<<" "<<mmin2<<" "<<mmin3<<endl;
		
		ofstream maxs((string("maxs_")+string(name)+string(".out")).c_str());
		maxs<<mmin1<<" "<<mmin2<<" "<<mmin3<<endl;
		maxs<<mmin11<<" "<<mmin21<<" "<<mmin31<<endl;
		maxs<<mmin12<<" "<<mmin22<<" "<<mmin32<<endl;
		maxs<<llmax<<endl;		
		//for(uint j=0;j<npoints;j++) 
        //for(uint i=0;i<npoints;i++){
		//	int binmax=limits4->GetBin(i+1,j+1);
		//	maxs<<"("<<i<<","<<j<<"):"<<limits4->GetBinContent(binmax)<<endl;
		//	}
    
	maxs.close();
		
    limits4->SetStats(0);
    limits4->GetXaxis()->SetTitle("log_{10}(tan\\beta)");
    limits4->GetYaxis()->SetTitle("M_{H+} (GeV)");
    
    Bmumu_Bsmumu->SetStats(0);
    Bmumu_Bsmumu->GetXaxis()->SetTitle("Br(B\\to\\mu\\mu)/10^{-10}");
    Bmumu_Bsmumu->GetYaxis()->SetTitle("Br(B_{s}\\to\\mu\\mu)/10^{-9}");
    
    
    limits_MR_MI->SetStats(0);
    limits_MR_MI->GetYaxis()->SetTitle("M_{I} (GeV)");
    limits_MR_MI->GetXaxis()->SetTitle("M_{R} (GeV)");
    
    
    
    limits_MR_McH->SetStats(0);
    limits_MR_McH->GetYaxis()->SetTitle("M_{H+} (GeV)");
    limits_MR_McH->GetXaxis()->SetTitle("M_{R} (GeV)");
    
    limits_MI_McH->SetStats(0);
    limits_MI_McH->GetYaxis()->SetTitle("M_{H+} (GeV)");
    limits_MI_McH->GetXaxis()->SetTitle("M_{I} (GeV)");
    
    limits_tb_MR->SetStats(0);
    limits_tb_MR->GetXaxis()->SetTitle("log_{10}(tan\\beta)");
    limits_tb_MR->GetYaxis()->SetTitle("M_{R} (GeV)");
    
    
    gStyle->SetOptTitle(0);
   
   Double_t contours[3];
   contours[0] = 0.003;
   contours[1] = 0.05;
   contours[2] = 0.32;
   
   
    double ma=0,me=.2, x0=1,y0=120;
     
    limits4->SetContour(3, contours);
    //limits4->GetYaxis()->SetLabelOffset(0.02);
    limits4->GetYaxis()->SetLabelSize(0.08);
    limits4->GetYaxis()->SetTitleSize(0.08);
    limits4->GetYaxis()->SetTitleOffset(1.2);
    limits4->GetYaxis()->SetLimits(1,999);
    //limits4->GetXaxis()->SetLabelOffset(0.02);
    limits4->GetXaxis()->SetLabelSize(0.08);
    limits4->GetXaxis()->SetTitleSize(0.08);
    limits4->GetXaxis()->SetTitleOffset(1.2);
    limits4->GetXaxis()->SetLimits(-2.99,2.99);
    
   
    
    
    gStyle->SetPaperSize(10.,10.);
    
    TLatex l;
    l.SetTextSize(0.08);
    string ss=qq[qup][gQ]+","+ll[lup][gL];
    
    TCanvas * c21=new TCanvas("c21","",int(800*(1+ma+me)),int(600*(1+ma+me)));
	c21->SetMargin(me,ma,me,ma);
	c21->SetGrid();
	//limits4->Draw("CONT Z LIST");
    limits4->Draw("CONT Z LIST");
    if(!mssm) l.DrawLatex(x0,y0,ss.c_str());
   
    c21->SaveAs((string("pdf_")+string(name)+string(".png")).c_str());
    
    delete c21;
    
    // Bmumu_Bsmumu->SetContour(3, contours);
    //limits4->GetYaxis()->SetLabelOffset(0.02);
    Bmumu_Bsmumu->GetYaxis()->SetLabelSize(0.08);
    Bmumu_Bsmumu->GetYaxis()->SetTitleSize(0.08);
    Bmumu_Bsmumu->GetYaxis()->SetTitleOffset(1.2);
    Bmumu_Bsmumu->GetYaxis()->SetLimits(0.01,4.99);
    //limits4->GetXaxis()->SetLabelOffset(0.02);
    Bmumu_Bsmumu->GetXaxis()->SetLabelSize(0.08);
    Bmumu_Bsmumu->GetXaxis()->SetTitleSize(0.08);
    Bmumu_Bsmumu->GetXaxis()->SetTitleOffset(1.2);
    Bmumu_Bsmumu->GetXaxis()->SetLimits(0.01,1.99);
    
    TCanvas * cB=new TCanvas("cB","",int(800*(1+ma+me)),int(600*(1+ma+me)));
	cB->SetMargin(me,ma,me,ma);
	cB->SetGrid();
	//limits4->Draw("CONT Z LIST");
    Bmumu_Bsmumu->Draw("COLZ");
    l.DrawLatex(x0,y0,ss.c_str());
   
    cB->SaveAs((string("Bmumu_Bsmumu_")+string(name)+string(".png")).c_str());
    
    delete cB;
    
    
    limits_MR_MI->SetContour(3, contours);
    //limits4->GetYaxis()->SetLabelOffset(0.02);
    limits_MR_MI->GetYaxis()->SetLabelSize(0.06);
    limits_MR_MI->GetYaxis()->SetTitleSize(0.06);
    limits_MR_MI->GetYaxis()->SetTitleOffset(1.1);
    limits_MR_MI->GetYaxis()->SetLimits(1,999);
    //limits4->GetXaxis()->SetLabelOffset(0.02);
    limits_MR_MI->GetXaxis()->SetLabelSize(0.06);
    limits_MR_MI->GetXaxis()->SetTitleSize(0.06);
    limits_MR_MI->GetXaxis()->SetTitleOffset(1.1);
    limits_MR_MI->GetXaxis()->SetLimits(1,999);
    
	TCanvas * c3=new TCanvas("c3","",800,600);
	c3->SetMargin(me,ma,me,ma);
	c3->SetGrid();
	
	limits_MR_MI->Draw("CONT LIST");

    c3->SaveAs((string("pdf_")+string(name)+string("_MRMI.png")).c_str());
	
	limits_MR_McH->SetContour(3, contours);
	//limits4->GetYaxis()->SetLabelOffset(0.02);
    limits_MR_McH->GetYaxis()->SetLabelSize(0.06);
    limits_MR_McH->GetYaxis()->SetTitleSize(0.06);
    limits_MR_McH->GetYaxis()->SetTitleOffset(1.1);
    limits_MR_McH->GetYaxis()->SetLimits(1,999);
    //limits4->GetXaxis()->SetLabelOffset(0.02);
    limits_MR_McH->GetXaxis()->SetLabelSize(0.06);
    limits_MR_McH->GetXaxis()->SetTitleSize(0.06);
    limits_MR_McH->GetXaxis()->SetTitleOffset(1.1);
    limits_MR_McH->GetXaxis()->SetLimits(1,999);
    
	TCanvas * c4=new TCanvas("c4","",800,600);
	limits_MR_McH->Draw("CONT LIST");	
    c4->SetMargin(me,ma,me,ma);
	c4->SetGrid();
    c4->SaveAs((string("pdf_")+string(name)+string("_MRMcH.png")).c_str());
	
	limits_MI_McH->SetContour(3, contours);
	//limits4->GetYaxis()->SetLabelOffset(0.02);
    limits_MI_McH->GetYaxis()->SetLabelSize(0.06);
    limits_MI_McH->GetYaxis()->SetTitleSize(0.06);
    limits_MI_McH->GetYaxis()->SetTitleOffset(1.1);
    limits_MI_McH->GetYaxis()->SetLimits(1,999);
    //limits4->GetXaxis()->SetLabelOffset(0.02);
    limits_MI_McH->GetXaxis()->SetLabelSize(0.06);
    limits_MI_McH->GetXaxis()->SetTitleSize(0.06);
    limits_MI_McH->GetXaxis()->SetTitleOffset(1.1);
    limits_MI_McH->GetXaxis()->SetLimits(1,999);
    
	TCanvas * c6=new TCanvas("c6","",800,600);
	limits_MI_McH->Draw("CONT LIST");	
    c6->SetMargin(me,ma,me,ma);
	c6->SetGrid();
    c6->SaveAs((string("pdf_")+string(name)+string("_MIMcH.png")).c_str());
	
	TCanvas * c5=new TCanvas("c5","",800,600);
	limits_tb_MR->Draw("colz");	

    c5->SaveAs((string("pdf_")+string(name)+string("_tbMR.png")).c_str());
	
	
	cout<<"gL "<<gL<<" gQ "<<gQ<<endl;
	cout<<"lup "<<lup<<" qup "<<qup<<endl;
	cout<<"llmax "<<llmax<<" gmax "<<gmax<<endl;


	delete m;
	delete limits4;
	delete Bmumu_Bsmumu;	
	delete limits_tb_MR;
	delete limits_tb_MI;
	delete limits_MR_MI;
	delete limits_MI_McH;
	delete limits_MR_McH;
	
    //mass.close();
	return 0;
	
}

