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
#include "TH2F.h"
#include "TVector.h"
#include "TCanvas.h"
#include "TMath.h"
#include <iostream>
#include <fstream>


#include <cln/cln.h>
#include <cln/float.h>

using namespace BGLmodels;

/**
* @brief A second implementation of the BGL model, for testing purposes
*/
class BGL2: public Model{
public:

BGL2(int genL=2,int genQ=2, int lup=0, int qup=0, int mssm=0): 
	   planck(6.58211928e-25),
	   GF("G_F"),
	   MZ("M_Z"),
	   MW("M_W"),
	   Mpip("Mpip",0.1396,"M_{\\pi^+}",domain::real),
	   Mpi0("Mpi0",0.1349766,"M_{\\pi^0}",domain::real),
	   MBp("MBp",5.279,"M_{B^+}",domain::real),
	   MB0("MB0",5.2795,"M_{B^0}",domain::real),
	   MBs0("MBs0",5.3663,"M_{B_s^0}",domain::real),
	   MKp("MKp",0.493677,"MKp",domain::real),
	   MK0("MK0",0.497614,"MK0",domain::real),
	   MDp("MDp",1.86957,"MDp",domain::real),
	   MD0("MD0",1.86480,"MD0",domain::real),
	   MDsp("MDsp",1.96845,"MDsp",domain::real),
	   MDs0("MDs0",0),
	   Fpi("Fpi",0.132,"Fpi",domain::real), 
	   FB("FB",0.189,"FB",domain::real),
	   FBs("FBs",0.225,"FBs",domain::real), 
	   FK("FK",0.159,"FK",domain::real),
	   FD("FD",0.208,"FD",domain::real),
	   FDs("FDs",0.248,"FDs",domain::real),
	   //alpha(7.297352e-3*4*M_PI),
	   cos2(pow(MW/MZ,2)),
	   g(sqrt(GF*8/sqrt(ex(2)))*MW),
	   //g(sqrt(4*Pi*alpha/(1-cos2))),
	   tanb("tg\\beta"),
	   cp("cp"),
	   McH("M_{H^+}"),
       MR("M_{R}"),
       MI("M_{I}"),
       Tparam("T_param"),
       Sparam("S_param"),
       QCD1("QCD_1"),
       QCD2("QCD_2"),
	   mixes(tanb,cp, genL,genQ, lup, qup, mssm),
	   mu("\\mu"),
	   BGLtype(4,0),
	   mmmax(1000),
	   stepsize(1e-2)
	   {
	alpha=pow(g,2)*(1-cos2)/(4*Pi);
	replacements.append(GF==1.166371e-5);
	replacements.append(MZ==M_MZ);
	replacements.append(MW==M_MW);
		
    mixes.appendtolst(replacements);
    
    replacements.append(Pi==M_PI);
    replacements.append(sqrt(ex(2))==sqrt(2));
        replacements.append(Pi==M_PI);
    replacements.append(sqrt(ex(2))==sqrt(2));
   
	Boson boson;
	
	realsymbol q3("q3");
	ex vq3=dirac_slash(q3,4);
	varidx jmu(mu,4,1);
	
	for(uint i=0;i<2;i++)
		for(uint j=0;j<3;j++)
			for(uint k=0;k<3;k++){
				conjtoabs.append(conjugate(mixes.V[i][j][k])==pow(abs(mixes.V[i][j][k]),2)/mixes.V[i][j][k]);
				}
	
	//W+ boson
	boson.mass=MW;
	boson.s=sVector;
	
	for(uint t=tLepton;t<=tQuark;t++) boson.C[t][iUp][iDown][hLeft]=mixes.V[t]*Matrix(g/sqrt(ex(2)));
	Boson wboson=boson;
	bosons.push_back(boson);
	boson.reset();
	
	//H+ boson
	boson.mass=McH;
	boson.s=sScalar;
	
	for(uint t=tLepton;t<=tQuark;t++)
	for(uint i=iUp;i<=iDown;i++) boson.C[t][iUp][iDown][i]=mixes.VN[t][i]*Matrix(g/MW/sqrt(ex(2)));
	Boson chiggs=boson;
	bosons.push_back(boson);
	boson.reset();
	
	for(int b=bosons.size()-1;b>=0;b--){
		boson.mass=bosons[b].mass;
		boson.s=bosons[b].s;
		if(boson.s==sVector)
			for(uint t=tLepton;t<=tQuark;t++)
			for(uint i=iUp;i<=iDown;i++)
			for(uint j=iUp;j<=iDown;j++)
			for(uint h=hLeft;h<=hRight;h++){
				boson.C[t][i][j][h]=bosons[b].C[t][j][i][h].conjugate();
				}	
		else for(uint t=tLepton;t<=tQuark;t++)
			for(uint i=iUp;i<=iDown;i++)
			for(uint j=iUp;j<=iDown;j++)
			for(uint h=hLeft;h<=hRight;h++){
				boson.C[t][i][j][hLeft]=bosons[b].C[t][j][i][hRight].conjugate();
				boson.C[t][i][j][hRight]=bosons[b].C[t][j][i][hLeft].conjugate();
				}			
		bosons.push_back(boson);
		boson.reset();
		}
	
	//(R+iI)/sqrt(2) boson
	boson.mass=MR;
	boson.s=sScalar;
	
	for(uint t=tLepton;t<=tQuark;t++){
			boson.C[t][iDown][iDown][hRight]=mixes.N[t][iDown]*Matrix(g/MW/ex(2));
			boson.C[t][iUp][iUp][hLeft]=mixes.N[t][iUp].conjugate()*Matrix(g/MW/ex(2));
			boson.C[t][iDown][iDown][hLeft]=mixes.N[t][iDown].conjugate()*Matrix(g/MW/ex(2));
			boson.C[t][iUp][iUp][hRight]=mixes.N[t][iUp]*Matrix(g/MW/ex(2));		
		}
	bosons.push_back(boson);
	boson.reset();
	
	//(R+iI)/sqrt(2) boson
	boson.mass=MI;
	boson.s=sScalar;
	
	for(uint t=tLepton;t<=tQuark;t++){
			boson.C[t][iDown][iDown][hRight]=mixes.N[t][iDown]*Matrix(I*g/MW/ex(2));
			boson.C[t][iUp][iUp][hLeft]=mixes.N[t][iUp].conjugate()*Matrix(I*g/MW/ex(2));
			boson.C[t][iDown][iDown][hLeft]=mixes.N[t][iDown].conjugate()*Matrix(-I*g/MW/ex(2));
			boson.C[t][iUp][iUp][hRight]=mixes.N[t][iUp]*Matrix(-I*g/MW/ex(2));		
		}
	bosons.push_back(boson);
	boson.reset();
	
	Fermion electron(tLepton,iDown,fElectron);
	Fermion electronR(tLepton,iDown,fElectron,cParticle,hRight);
	
	Fermion muon(tLepton,iDown,fMuon);
	Fermion muonR(tLepton,iDown,fMuon,cParticle,hRight);
	
	Fermion tau(tLepton,iDown,fTau);
	Fermion tauR(tLepton,iDown,fTau,cParticle,hRight);
	Fermion neutrino(tLepton,iUp);
	Fermion neutrinotau(tLepton,iUp,fTau);
	Fermion neutrinomuon(tLepton,iUp,fMuon);
	Fermion neutrinoe(tLepton,iUp,fElectron);
	
	Fermion up(tQuark,iUp,fElectron);
	Fermion down(tQuark,iDown,fElectron);
	Fermion bottom(tQuark,iDown,fTau);
	Fermion strange(tQuark,iDown,fMuon);
	Fermion charm(tQuark,iUp,fMuon);
	Fermion top(tQuark,iUp,fTau);
	
	Meson Pi0d(down,down,Mpi0,Fpi);
	Meson Pi0u(down,down,Mpi0,Fpi);
	Meson Pip(up,down,Mpip,Fpi);
	Meson Pim(down,up,Mpip,Fpi);
	
	Meson K0(down,strange,MK0,FK);
	Meson Kp(up,strange,MKp,FK);
	
	Meson D0(charm,up,MD0,FD);
	Meson Dp(charm,down,MDp,FD);
	Meson Dsp(charm,strange,MDsp,FDs);
	
	Meson B0(down,bottom,MB0,FB);
	Meson Bp(up,bottom,MBp,FB);
	Meson Bs0(strange,bottom,MBs0,FBs);

	lst sb;
	//sb.append(mixes.M[tQuark][iUp][0][0]==0);
	sb.append(pow(abs(mixes.V[0][2][2]),2)==1-pow(abs(mixes.V[0][1][2]),2)-pow(abs(mixes.V[0][0][2]),2));
	sb.append(pow(abs(mixes.V[0][2][1]),2)==1-pow(abs(mixes.V[0][1][1]),2)-pow(abs(mixes.V[0][0][1]),2));
	
	//cout<<pow(sqrt(2)/8*pow(g/MW,2),2)<<endl;
	//cout<<pow(1.166,2)<<endl;
	 double fK=0.156;
    ex KKbar=ex(std::pow(fK,2))*mesonmixing(MK0,strange,down);
    ex eK=0.94*imag_part(KKbar)/3.5e-15/sqrt(2);
    cout<<"KKbar "<<KKbar<<endl;
    KKbar=expand(KKbar.subs(replacements).subs(lst(abs(wild()*pow(MR,-2))==abs(wild())*pow(MR,-2))).subs(lst(log(wild()*pow(MR,-2))==log(wild())-2*log(MR))));
	KKbar=expand(KKbar.evalf());
	ex aKKbar=sqrt(KKbar.real_part()*KKbar.real_part()+KKbar.imag_part()*KKbar.imag_part());   
   
    eK=0.94*imag_part(KKbar)/3.5e-15/sqrt(2);
    eK=eK.subs(replacements).real_part();
	eK=collect_common_factors(expand(eK.evalf()));
	cout<<"eK"<<eK<<endl;
    //add("a_eK",abs(eK),new limitedobs(2.2e-3));
    epsilonK=new calcuex(new limitedobs(2*0.011e-3),abs(eK));
    
    
}

~BGL2(){epsilonK->~calcuex();}

parameters generateparameters(int max=0) const{
	parameters p;
	//x=log_10(tanb)
	p.push_back(freeparameter(-3,3,r,stepsize));
	//y=log_10(McH)
	if(max==1) p.push_back(freeparameter(10,10000,r,stepsize));
	else p.push_back(freeparameter(10,mmmax,r,stepsize));
	//log_10(massR)
	p.push_back(freeparameter(-200,200,r,stepsize));
	//log_10(massI)
	p.push_back(freeparameter(-50,50,r,stepsize));
	
	return p;
}


parameters getlist(const parameters & p) const{
	//cout<<aux<<endl;
	//double c2=(1+sqrt(1-4*sqrt(ex_to<numeric>(mudecay.subs(lst(tanb==exp(p[0].value),McH==p[1].value))).to_double())))/2;
	
	double x=pow(10.0,p[0].value);
	//double y=pow(10.0,p[1].value);
	//double z=pow(10.0,p[2].value);
	//double w=pow(10.0,p[3].value);
	
	double y=p[1].value;
	double z=y+p[2].value;
	double w=z+p[3].value;
	
	parameters pp(p);
	pp[0].value=x;
	pp[2].value+=pp[1].value;
	pp[3].value+=pp[2].value;
	pp.values=vector<double>();
	for(uint i=0; i<4; i++) pp.values.push_back(pp[i].value);
	lst &l=pp.p;
	l=lst(tanb==x,McH==y,MR==z,MI==w);
    l.append(QCD1==inter1.Eval(y));
	l.append(QCD2==inter2.Eval(y));
	
	for(uint i=0;i<3;i++){
		l.append(Mu[i]==Mu_[i].Eval(log(y)));
		l.append(Md[i]==Md_[i].Eval(log(y)));
	}
	return pp;
}

ex mesonmixing(ex mesonmass, const Fermion& f1, const Fermion& f2) const{
	
	ex ret=0;
	
		ex v1=0, v2=0;
		ex mq1=mixes.mass(f1),mq2=mixes.mass(f2);
		ex m2q1=mq1*mq1, m2q2=mq2*mq2;
	
	for(uint i=0;i<bosons.size();i++)
		if(bosons[i].s==0){
			ex a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1));
			v1=v1+pow(a/bosons[i].mass,2);
			
			ex b=(bosons[i].couplingdaggerR(f2,f1)+bosons[i].couplingdaggerL(f2,f1));
			v2=v2+pow(b/bosons[i].mass,2);	
		}
	
	ex fc=mesonmass/(mixes.massnum(f1)+mixes.massnum(f2));
	fc=pow(fc,2);
	
	ret=2*(-v1*(1+11*fc)+v2*(1+fc))*mesonmass/96;
	
	return collect_common_factors(ret.subs(conjtoabs));
	//return expand(ret.subs(lst(exp(-I*wild())==1/exp(I*wild()),sin(wild())==sqrt(1-pow(cos(wild()),2)))));
}

double bsgammawidth(double tanb_,double McH_,double MR_,double MI_, int option=0){
	parameters p=generateparameters();
	p[0].value=pow(10.0,tanb_);
	p[1].value=McH_;
	p[2].value=MR_;
	p[3].value=MI_;
	calcubtosgamma2 cal(mixes);
	
	return cal.width(p,option);
}

double epsK(double tanb_,double McH_,double MR_,double MI_, int option=0){
	parameters p=generateparameters();
	p[0].value=pow(10.0,tanb_);
	p[1].value=McH_;
	p[2].value=MR_;
	p[3].value=MI_;
    p.p=lst(tanb==p[0].value,McH==p[1].value,MR==p[2].value,MI==p[3].value);

	 
	 return epsilonK->error(p);
}

const double planck;
const possymbol GF, MZ, MW, Mh;
const constant Mpip, Mpi0, MBp,MB0,MBs0, MKp,MK0,MDp,MD0,MDsp,MDs0;
const constant Fpi, FB,FBs, FK,FD,FDs;
ex cos2, g, alpha;
const possymbol tanb, cp, McH, MR, MI, rho;
const realsymbol Tparam, Sparam, QCD1, QCD2;
possymbol Mu[3],Md[3];
vector< Boson > bosons;

lst replacements;
ex Btaunu;
ex BR_Htotaunu;
ex BR_toptoHq;
ex BtotaunuR;
ex BtoDtaunuR;
ex BtoD2taunuR;

const Mixes mixes;
lst conjtoabs;
realsymbol mu;

int iBtaunu, iBDtaunu, iBD2taunu;
vector<int> BGLtype;
ROOT::Math::Interpolator inter1, inter2;  
ROOT::Math::Interpolator Mu_[3],Md_[3];
double mmmax, stepsize;

calcuex * epsilonK;

};


/**
* @brief the main function takes the arguments inputfile gL gQ lup qup which specify the file containing the simulation results for a BGL model and draws the plots for that model 
*/
int main(int argc, char* argv[]){
    // Check the number of parameters
    
    if(argc<6){
        std::cerr<<"Usage: "<<argv[0]<<" inputfile gL gQ lup qup"<<std::endl;
        return 1;}
    CD cmu=conj(Vud[1][0])*Vud[1][1];
    CD umu=conj(Vud[0][0])*Vud[0][1];
    
    
   cout<<"RATIO "<<cmu*cmu<<endl;
   cout<<umu*umu<<endl;

   int gL=atoi(argv[2]);
   int gQ=atoi(argv[3]);
   int lup=atoi(argv[4]);
   int qup=atoi(argv[5]);
   char name[5]="0000";
	
   name[0]+=gL;
   name[1]+=gQ;
   name[2]+=lup;
   name[3]+=qup;
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
    double initBmumu=0, finalBmumu=3;
    double initBsmumu=0, finalBsmumu=6;
    
    double llmax=-1000,McHmax=1000,MRmax=1000,MImax=1000,tbmax=1;
	
	TFile *f=new TFile(argv[1],"read");
	if(!f->IsOpen()) cout<<"NOFILE"<<endl;
	//f->ShowStreamerInfo();

	TH2F *limits4,*Bmumu_Bsmumu,*limits_tb_MR,*limits_tb_MI;
	TH2F *limits_MR_MI,*limits_MR_McH,*limits_MI_McH;
	
	f->GetObject("limits4;1",limits4);
	f->GetObject("Bmumu_Bsmumu;1",Bmumu_Bsmumu);
	f->GetObject("limits_tb_MR;1",limits_tb_MR);
	f->GetObject("limits_tb_MI;1",limits_tb_MI);
	f->GetObject("limits_MR_MI;1",limits_MR_MI);
	f->GetObject("limits_MR_McH;1",limits_MR_McH);
	f->GetObject("limits_MI_McH;1",limits_MI_McH);
	
	TVectorD* vllmax=NULL;
	
	f->GetObject("vllmax;1",vllmax);
	if(!vllmax) cout<<"ERROR"<<endl;
	llmax=(*vllmax)[0];
	//tbmax=(*vllmax)[1];
	//McHmax=(*vllmax)[2];
	//MRmax=McHmax+(*vllmax)[3];
	//MImax=MRmax+(*vllmax)[4];
	cout<<llmax<<" "<<tbmax<<" "<<McHmax<<" "<<MRmax<<" "<<MImax<<" "<<endl;	
	
  /*BGL2* m=new BGL2(gL,gQ,lup,qup);
  double sm_(m->bsgammawidth(tbmax,McHmax,MRmax,MImax,5));
  double charged_(m->bsgammawidth(tbmax,McHmax,MRmax,MImax,1));
  double neutral_(m->bsgammawidth(tbmax,McHmax,MRmax,MImax,2));
  double neutralR_(m->bsgammawidth(tbmax,McHmax,MRmax,MImax,3));
  double neutralI_(m->bsgammawidth(tbmax,McHmax,MRmax,MImax,4));
  double eK_(m->epsK(tbmax,McHmax,MRmax,MImax));
  
  double all_(m->bsgammawidth(tbmax,McHmax,MRmax,MImax,0));
  */	
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
			//int nn=4;
			//int ii=(i/nn)*nn, jj=(j/nn)*nn;
			//for(int iii=ii;iii<ii+n;++iii)
			//for(int iii=ii;iii<ii+n;++iii)
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
		maxs<<llmax<<" "<<tbmax<<" "<<McHmax<<" "<<MRmax<<" "<<MImax<<" "<<endl;
		//maxs<<eK_<<endl;
		//maxs<<sm_<<" "<<charged_<<" "<<neutral_<<" "<<neutralR_<<" "<<neutralI_<<" "<<all_<<endl;
				
		//for(uint j=0;j<npoints;j++) 
        //for(uint i=0;i<npoints;i++){
		//	int binmax=limits4->GetBin(i+1,j+1);
		//	maxs<<"("<<i<<","<<j<<"):"<<limits4->GetBinContent(binmax)<<endl;
		//	}
    
	maxs.close();
	
	   double ma=0,me=.2, x0=1,y0=120;
   gStyle->SetOptTitle(0);
   gStyle->SetPaperSize(10.,10.);
   TCanvas * c21=new TCanvas("c21","",int(800*(1+ma+me)),int(600*(1+ma+me)));
	c21->SetMargin(me,ma,me,ma);
	c21->SetGrid();
	
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
    
    
  
   Double_t contours[3];
   contours[0] = 0.003;
   contours[1] = 0.05;
   contours[2] = 0.32;
   
   
 
     
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
    
    
    
    TLatex l;
    l.SetTextSize(0.08);
    string ss=qq[qup][gQ]+","+ll[lup][gL];
    
    
	limits4->Draw("CONT Z LIST");
    //limits4->Draw("CONT LIST");
    //limits4->Draw("colz");
    
    l.DrawLatex(x0,y0,ss.c_str());
   
    c21->SaveAs((string("pdf_")+string(name)+string(".png")).c_str());
    
    delete c21;
    //Bmumu_Bsmumu->SetBit(TH1::kCanRebin);
    Bmumu_Bsmumu->Rebin2D(2,2);
     
    //Bmumu_Bsmumu->SetContour(3, contours);
    //limits4->GetYaxis()->SetLabelOffset(0.02);
    Bmumu_Bsmumu->GetYaxis()->SetLabelSize(0.08);
   Bmumu_Bsmumu->GetYaxis()->SetTitleSize(0.08);
    Bmumu_Bsmumu->GetYaxis()->SetTitleOffset(0.8);
    Bmumu_Bsmumu->GetYaxis()->SetLimits(0.01,4.99);
    //Bmumu_Bsmumu->GetYaxis()->SetRangeUser(0.01, 3.49);
   // Bmumu_Bsmumu->GetYaxis()->SetLimits(0.01,3.49);
    //limits4->GetXaxis()->SetLabelOffset(0.02);
    Bmumu_Bsmumu->GetYaxis()->SetNdivisions(5, kTRUE);
    Bmumu_Bsmumu->GetXaxis()->SetNdivisions(5, kTRUE);
    
    Bmumu_Bsmumu->GetXaxis()->SetLabelSize(0.08);
    Bmumu_Bsmumu->GetXaxis()->SetTitleSize(0.08);
   Bmumu_Bsmumu->GetXaxis()->SetTitleOffset(1.2);
    Bmumu_Bsmumu->GetXaxis()->SetLimits(0.01,1.99);
    //Bmumu_Bsmumu->GetXaxis()->SetRangeUser(0., 2);
     
     TCanvas * cB=new TCanvas("cB","",int(800*(1+ma+me)),int(600*(1+ma+me)));
	cB->SetMargin(.14,ma,me,ma);
	cB->SetGrid();
	//limits4->Draw("CONT Z LIST");
    Bmumu_Bsmumu->Draw("COLZ");
    
    l.DrawLatex(1.5,1,ss.c_str());
   
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
	
    //delete m;
    //mass.close();
    f->Close();
	delete f;
	return 0;
	
}

