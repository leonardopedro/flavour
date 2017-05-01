#ifndef BGL_H
#define BGL_H

#include "widthcalc.h"

#include "TH2F.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include <iostream>

#include "Math/Polynomial.h"
#include "Math/Interpolator.h"
#include "Formulas.h"
 

namespace BGLmodels{

/**
* @brief Gauge boson
*/
class Boson {
public:
	
	Boson(): C(Matrix(),2,2,2,2){}
    
	ex couplingL(const Fermion& f2,const Fermion& f1) const { 
		bool quiralfilter=0;
		if(f1.type!=f2.type) return 0;
		if(s==sScalar){
			if(f1.helicity!=hRight && f2.helicity!=hLeft) quiralfilter=1;
			}
		else if(s==sVector){
			if(f1.helicity!=hRight && f2.helicity!=hRight) quiralfilter=1;
			}
		
		if(quiralfilter) return C[f2.type][f2.isospin][f1.isospin][hLeft][f2.flavour][f1.flavour];
		return 0;
		}
	ex couplingR(const Fermion& f2,const Fermion& f1) const { 
		bool quiralfilter=0;
		if(f1.type!=f2.type) return 0;
		if(s==sScalar){
			if(f2.helicity!=hRight && f1.helicity!=hLeft) quiralfilter=1;
			}
		else if(s==sVector){
			if(f1.helicity!=hLeft && f2.helicity!=hLeft) quiralfilter=1;
			}
		if(quiralfilter) return C[f2.type][f2.isospin][f1.isospin][hRight][f2.flavour][f1.flavour];
		return 0;
		}
	ex couplingdaggerL(const Fermion& f2,const Fermion& f1) const {
		if(s==sScalar) return couplingR(f1,f2).conjugate();
		return couplingL(f1,f2).conjugate();
		}
	ex couplingdaggerR(const Fermion& f2,const Fermion& f1) const {
		if(s==sScalar) return couplingL(f1,f2).conjugate();
		return couplingR(f1,f2).conjugate();
		}
	ex coupling(const Fermion& f2,const Fermion& f1, ex mu){
		if(s==sScalar) return couplingL(f2,f1)*dirac_gammaL()+couplingR(f2,f1)*dirac_gammaR();
		else return couplingL(f2,f1)*dirac_gammaL()+couplingR(f2,f1)*dirac_gammaR();
		}
	void reset(){
		C=multivector<Matrix,4>(Matrix(),2,2,2,2);
		}
    BSpin s;
    ex mass;
    multivector<Matrix,4> C;
};


	
//Boson::Type Boson::dagger[2][2]={Boson::scalarright,Boson::scalarleft,Boson::vectorleft,Boson::vectorright};

/**
* @brief Implementation of the BGL model
*/
class BGL: public Model{
public:

BGL(int genL=2,int genQ=2, int lup=0, int qup=0, int mssm=0): 
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
	   mixes(tanb,cp, genL,genQ, lup, qup, mssm),
	   mu("\\mu"),
	   BGLtype(4,0),
	   mmmax(1000),
	   stepsize(1e-2)
	   //muwidth(planck/2.197034e-6)
	   {
	alpha=pow(g,2)*(1-cos2)/(4*Pi);
	replacements.append(GF==1.166371e-5);
	replacements.append(MZ==M_MZ);
	replacements.append(MW==M_MW);
		
    mixes.appendtolst(replacements);
    
    replacements.append(Pi==M_PI);
    replacements.append(sqrt(ex(2))==sqrt(2));
    
	//cout<<pow(sqrt(2)/8*pow(g/MW,2),2)<<endl;
	//cout<<pow(1.166,2)<<endl;
	
	Boson boson;
	
	realsymbol q3("q3");
	ex vq3=dirac_slash(q3,4);
	varidx jmu(mu,4,1);
	
	for(uint i=0;i<2;i++)
		for(uint j=0;j<3;j++)
			for(uint k=0;k<3;k++){
				conjtoabs.append(conjugate(mixes.V[i][j][k])==pow(abs(mixes.V[i][j][k]),2)/mixes.V[i][j][k]);
				}
	/*
	//Gamma boson
	boson.mass=0;
	boson.s=Boson::vector;
	
	boson.coupsL[0][0]=Matrix(g*sqrt(1-cos2)*0);
	boson.coupsL[1][1]=Matrix(g*sqrt(1-cos2)*(-1));
	boson.coupsL[2][2]=Matrix(g*sqrt(1-cos2)*ex(2)/3);
	boson.coupsL[3][3]=Matrix(g*sqrt(1-cos2)*ex(-1)/3);
	
	boson.coupsR[0][0]=Matrix(g*sqrt(1-cos2)*0);
	boson.coupsR[1][1]=Matrix(g*sqrt(1-cos2)*(-1));
	boson.coupsR[2][2]=Matrix(g*sqrt(1-cos2)*ex(2)/3);
	boson.coupsR[3][3]=Matrix(g*sqrt(1-cos2)*ex(-1)/3);
	
	bosons.push_back(boson);
	boson.reset();
	*/
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
	
	//cout<<"Btaunu "<<collect_common_factors(expand(Btaunu.subs(sb).subs(conjtoabs)))<<endl;
	
	cout<<latex;
	
	ex mutoenunu=decaywidth(muon,neutrino,electron,neutrino);
	
	//cout<<"mutoenunu "<<mutoenunu<<endl;
	//add("mutoenunu",decaywidth(muon,neutrino,electron,neutrino),new gaussobs(planck/2.197034e-6,0.03));
	
	add("muRtoeRnunu",gRR2(muon,electron),new limitedobs(std::pow(0.035,2),0.95));
	
	//add("tautoenunu",decaywidth(tau,neutrino,electron,neutrino),new gaussobs(planck/290.6e-15*0.1782,0.03));
	add("tauRtoeRnunu",gRR2(tau,electron),new limitedobs(std::pow(0.7,2),0.95));
	
	//add("tautomununu",decaywidth(tau,neutrino,muon,neutrino),new gaussobs(planck/290.6e-15*0.1739,0.03));
	add("tauRtomuRnunu",gRR2(tau,muon),new limitedobs(std::pow(0.72,2),0.95));
	
	add("tautomu_tautoe",tautomu_tautoe(),new gaussobs(1.0018,0.0014/1.0018)); //PROBLEM!!!
	cout<<"tautomu_tautoe: "<<1/1.0018<< "ERROR: "<<0.0014/1.0018<<endl;
	cout<<"ratio1 "<<tautomu_tautoe().subs(replacements)<<endl;
	cout<<"ratio2 "<<(decaywidth(tau,neutrino,muon,neutrino,sVector)/decaywidth(tau,neutrino,electron,neutrino,sVector)).subs(replacements)<<endl;
	//muto3e
	ex mu3e=decaywidth(muon,electron,electron,electron);
	add("muto3e", mu3e,new limitedobs(planck/2.197034e-6*1e-12));
	cout<<"mu3e "<<decaywidthtest2(muon)<<endl;
	
	//tauto3e
	add("tauto3e", decaywidth(tau,electron,electron,electron),new limitedobs(planck/290.6e-15*2.7e-8));
	//tauto2e1mu+
	add("tauto2e1mup", decaywidth(tau,electron,electron,muon), new limitedobs(planck/290.6e-15*1.5e-8));
	//tauto2e1mu-
	add("tauto2e1mu", decaywidth(tau,electron,muon,electron), new limitedobs(planck/290.6e-15*1.8e-8));
	//tauto2mu1e+
	
	add("tauto2mu1ep", decaywidth(tau,muon,muon,electron), new limitedobs(planck/290.6e-15*1.7e-8));
	cout<<"tauto2mu1ep "<<decaywidthtest2(tau)<<endl;
	//tauto2mu1e-
	add("tauto2mu1ep", decaywidth(tau,muon,electron,muon), new limitedobs(planck/290.6e-15*2.7e-8));
	cout<<"tauto2mu1e "<<decaywidthtest2(tau)<<endl;
	
	//tauto3mu
	add("tauto3mu", decaywidth(tau,muon,muon,muon),new limitedobs(planck/290.6e-15*2.1e-8));
	//cout<<"tauto3mu "<<collect_common_factors(expand(decaywidth(tau,muon,muon,muon)))<<endl;
	

	ex piratio=1.2352e-4/(mesondw(Pip,neutrino,electron,sVector)/mesondw(Pip,neutrino,muon,sVector));
	ex picorrection=piratio.subs(replacements);
	ex pierror=picorrection*0.0001/1.2352;
	cout<<"PiRatio "<<picorrection-1<<" +/- "<<pierror<<endl;
	piratio*=mesondw(Pip,neutrino,electron)/mesondw(Pip,neutrino,muon);
	add("piontoenu_munu",piratio,new gaussobs(1.230e-4,0.003)); //PROBLEM!!!
	cout<<"piontoenu_munu: "<<1.2352e-4/1.230e-4<<" ERROR: "<<0.003<<endl;
	
	
	add("tautopinu_pitomunu",(1+0.16e-2)*fermiontomeson(tau,neutrino,Pip)/mesondw(Pip,neutrino,muon),new gaussobs((10.83e-2/290.6e-15/(0.9998770/2.6033e-8)),0.06/10.83));
	cout<<"tautopinu/pitomunu: "<<(1+0.16e-2)*(fermiontomeson(tau,neutrino,Pip,sVector)/mesondw(Pip,neutrino,muon,sVector)).subs(replacements)/((10.83e-2/290.6e-15/(0.9998770/2.6033e-8)))<<" ERROR: "<<0.06/10.83<<endl;
	cout<<"tautopinu: "<<fermiontomesontest(tau,neutrino,Pip)<<endl;
	cout<<"tautopinu_pitomunu: "<<10.83e-2/290.6e-15/(0.9998770/2.6033e-8)<<" +/- "<<0.06e-2/290.6e-15/(0.9998770/2.6033e-8)<<endl;
	
	add("tautoKnu_Ktomunu",(1+0.9e-2)*fermiontomeson(tau,neutrino,Kp)/mesondw(Kp,neutrino,muon),new gaussobs((7e-3/290.6e-15)/(0.6355/1.238e-8),0.1/7));
    cout<<"tautoKnu/Ktomunu: "<<(1+0.9e-2)*(fermiontomeson(tau,neutrino,Kp,sVector)/mesondw(Kp,neutrino,muon,sVector)).subs(replacements)/((7e-3/290.6e-15)/(0.6355/1.238e-8))<<" ERROR: "<<0.1/7<<endl;
	cout<<"tautoKnu/Ktomunu: "<<(7e-3/290.6e-15)/(0.6355/1.238e-8)<<" +/- "<<(0.1e-3/290.6e-15)/(0.6355/1.238e-8)<<endl;
	
    //ex pi0toemu=(mesondecaywidth(Mpi0,down,down,electron,muon)+mesondecaywidth(Mpi0,up,up,electron,muon)+mesondecaywidth(Mpi0,down,down,muon,electron)+mesondecaywidth(Mpi0,up,up,muon,electron))/2;
    ex pi0toemu=(mesondw(Pi0d,electron,muon)+mesondw(Pi0d,muon,electron)+mesondw(Pi0u,electron,muon)+mesondw(Pi0u,muon,electron))/2;
    add("pi0toemu",pi0toemu,new limitedobs(3.6e-10*planck/8.52e-17));

	ex Kratio=2.477e-5/(mesondw(Kp,neutrino,electron,sVector)/mesondw(Kp,neutrino,muon,sVector));
	ex Kcorrection=Kratio.subs(replacements);
	ex Kerror=Kcorrection*0.001/2.477;
	cout<<"KRatio "<<Kcorrection-1<<" +/- "<<Kerror<<endl;
	Kratio*=mesondw(Kp,neutrino,electron)/mesondw(Kp,neutrino,muon);
	add("Ktoenu_munu",Kratio,new gaussobs(2.488e-5,0.005));
	cout<<"Ktoenu_munu: "<<2.477e-5/2.488e-5<<" ERROR: "<<0.005<<endl;
	
    ex k0Ltoemu=mesondw(K0,electron,muon)+mesondw(K0,muon,electron);
	add("K0Ltoemu",k0Ltoemu,new limitedobs((4.7e-12*planck/5.116e-8)));
    add("K0Ltoee",mesondw(K0,electron,electron),new limitedobs((9e-12*planck/5.116e-8)));
    
    //add("K0Ltomumu",mesondw(K0,muon,muon),new limitedobs((6.84e-9*planck/5.116e-8)));
  
	add("Dtoenu",mesondw(Dp,neutrino,electron),new limitedobs(8.8e-6*planck/1040e-15));
	add("Dtomunu",mesondw(Dp,neutrino,muon),new gaussobs(3.82e-4*planck/1040e-15,0.1)); //PROBLEM!!!
	cout<<"Dtomunu: "<<mesondw(Dp,neutrino,muon,sVector).subs(replacements)/(3.82e-4*planck/1040e-15)<<" ERROR: "<<0.1<<endl;
	
	add("Dtotaunu",mesondw(Dp,neutrino,tau),new limitedobs(1.2e-3*planck/1040e-15));  //PROBLEM!!!
	cout<<"Dtotaunu: "<<mesondw(Dp,neutrino,tau,sVector).subs(replacements)/(1.2e-3*planck/1040e-15)<<" LIMIT"<<endl;
	
	//D0 2.6e-7/410.1e-15
	
	ex D0toemu=mesondw(D0,electron,muon)+mesondw(D0,muon,electron);
	add("D0toemu",D0toemu,new limitedobs((2.6e-7*planck/410.1e-15)));
    ex D0toee=mesondw(D0,electron,electron);
	add("D0toee",D0toee,new limitedobs((7.9e-8*planck/410.1e-15)));
    ex D0tomumu=mesondw(D0,muon,muon);
	add("D0tomumu",D0tomumu,new limitedobs((1.4e-7*planck/410.1e-15)));
  
    //ex Dstomunu=mesondecaywidth(MDsp,strange,charm,muon,neutrino); //500e-15
	add("Dstomunu",mesondw(Dsp,neutrino,muon),new gaussobs(5.9e-3*planck/500e-15,0.33/5.9)); //PROBLEM!!!
	cout<<"Dstomunu: "<<mesondw(Dsp,neutrino,muon,sVector).subs(replacements)/(5.9e-3*planck/500e-15)<<" ERROR: "<<0.33/5.9<<endl;
	
	add("Dstoenu",mesondw(Dsp,neutrino,electron),new limitedobs(1.2e-4*planck/500e-15));
	add("Dstotaunu",mesondw(Dsp,neutrino,tau),new gaussobs(5.43e-2*planck/500e-15,0.31/5.43)); //PROBLEM!!!
	cout<<"Dstotaunu: "<<mesondw(Dsp,neutrino,tau,sVector).subs(replacements)/(5.43e-2*planck/500e-15)<<" ERROR: "<<0.31/5.43<<endl;
	
	add("Btomunu",mesondw(Bp,neutrino,muon),new limitedobs(9.8e-7*planck/1.641e-12));
	add("Btoenu",mesondw(Bp,neutrino,electron),new limitedobs(1e-6*planck/1.641e-12));
	
	add("Btotaunu",mesondw(Bp,neutrino,tau),new gaussobs(1.15e-4*planck/1.641e-12,0.23/1.15));
	//add("Btotaunu",mesondw(Bp,neutrino,tau),new gaussobs(0.79e-4*planck/1.641e-12,0.23/1.15));
	
	//calcuBmumu calcutest(mixes,Bs0,muon,muon,2.9e-9*planck/1.516e-12,0.7e-9*planck/1.516e-12,"Bs_to_mumu");
	//calcuBmumu calcutest2(mixes,B0,muon,muon,3.6e-10*planck/1.519e-12,1.6e-10*planck/1.516e-19,"B_to_mumu");
	//cout<<"TESTE "<<endl;
	//double ps[4]={1,1e16,1e16,1e16};
	//double resteste=0,resteste2=0;
	//int nt=4,mt=1;
	//calcutest.fp(&nt,ps,&mt,&resteste);
	//calcutest2.fp(&nt,ps,&mt,&resteste2);
	//cout<<"TESTE "<<resteste/(2.9e-9*planck/1.516e-12)<<" "<<resteste2/(3.6e-10*planck/1.519e-12)<<endl;
	//ex B0tomumu=mesondw(B0,muon,muon);
	//cout<<"B0tomumu "<<collect_common_factors(B0tomumu)<<endl;
	//1.65e-4
	//add("B0tomumu",B0tomumu,new limitedobs((8e-10*planck/1.519e-12)));

	push_back(prediction(new calcuBmumu(mixes,B0,muon,muon,new limitedobs(6.3e-10*planck/1.519e-12),"B_to_mumu")));
	//push_back(prediction(new calcuBmumu(mixes,B0,muon,muon,new gauss2obs(3.6e-10*planck/1.519e-12,1.6e-10*planck/1.519e-12),"B_to_mumu")));
	push_back(prediction(new calcuBmumu(mixes,Bs0,muon,muon,new gauss2obs(2.9e-9*planck/1.516e-12,0.7e-9*planck/1.516e-12),"Bs_to_mumu")));
    push_back(prediction(new calcuBmumu(mixes,K0,muon,muon,new limitedobs(2.3e-9*planck/5.116e-8),"K0L_to_mumu")));
    
    cBmumu=new calcuBmumu(mixes,B0,muon,muon,new limitedobs(6.3e-10*planck/1.519e-12),"B_to_mumu");
    cBsmumu=new calcuBmumu(mixes,Bs0,muon,muon,new gauss2obs(2.9e-9*planck/1.516e-12,0.7e-9*planck/1.516e-12),"Bs_to_mumu");
    
    ex B0toetau=mesondw(B0,electron,tau)+mesondw(B0,tau,electron);
	add("B0toetau",B0toetau,new limitedobs((2.8e-5*planck/1.519e-12)));
	ex B0tomutau=mesondw(B0,muon,tau)+mesondw(B0,tau,muon);
	add("B0tomutau",B0tomutau,new limitedobs((2.2e-5*planck/1.519e-12)));
	ex B0toee=mesondw(B0,electron,electron);
	add("B0toee",B0toee,new limitedobs((8.3e-8*planck/1.519e-12)));
    ex B0totautau=mesondw(B0,tau,tau);
	add("B0totautau",B0totautau,new limitedobs((4.1e-3*planck/1.519e-12)));
    
	//B0s m=5.3663, life=1.472e-12 emu=2e-7, ee=2.8e-7 mumu=4.2e-8
	ex Bs0toemu=mesondw(Bs0,electron,muon)+mesondw(Bs0,muon,electron);
	add("Bs0toemu",Bs0toemu,new limitedobs((2e-7*planck/1.516e-12)));
    ex Bs0toee=mesondw(Bs0,electron,electron);
	add("Bs0toee",Bs0toee,new limitedobs((2.8e-7*planck/1.516e-12)));
//    ex Bs0tomumu=mesondw(Bs0,muon,muon);
//	add("Bs0tomumu",Bs0tomumu,new limitedobs((3.2e-9*planck/1.516e-12)));
    
   // add("chargedHiggs",pow(McH,-2),new limitedobs(std::pow(80.0,-2),0.9));
    
    cout<<"Bs0tomumu: "<<mesondwtest(Bs0,muon,muon)<<endl;
    //add("chargedHiggs",1/McH,new limitedobs(1/80.0,0));
    
    /*
    Matrix llgamma2loop=Matrix(sqrt(ex(2))*mixes.N[tQuark][iUp][2][2]*mixes.M[tQuark][iUp][fTau][fTau]*pow(1/McH*log(mixes.M[tQuark][iUp][fTau][fTau]/McH),2))*mixes.N[tLepton][iDown];
    for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++)
			if(j<i) llgamma2loop[i][j]=ex(3)*pow(g*g*(1-cos2)/4/Pi/Pi,3)*llgamma2loop[j][i]*llgamma2loop[j][i].conjugate()/pow(mixes.M[tLepton][iDown][i][i],2);
			else llgamma2loop[i][j]=0;
	add("mutoegamma",llgamma2loop[1][0],new limitedobs(1.2e-11));
	add("tautoegamma",llgamma2loop[2][0],new limitedobs(3.3e-8));
	add("tautomugamma",llgamma2loop[2][1],new limitedobs(4.4e-8));
	//add("mutoegamma",llgamma2loop[1][0],new limitedobs(planck/2.197034e-6*1.2e-11));
	//add("tautoegamma",llgamma2loop[2][0],new limitedobs(planck/290.6e-15*3.3e-8));
	//add("tautomugamma",llgamma2loop[2][1],new limitedobs(planck/290.6e-15*4.4e-8));
    */
    
    Matrix llgammaCH;
	//Matrix llgammaCH=Matrix(g*g/(16*Pi*4*Pi*MW*MW*McH*McH*12))*mixes.M[tLepton][iDown]*mixes.VN[tLepton][1].conjugate()*mixes.VN[tLepton][1];
	Matrix llgammaH0M,llgammaH0E;
	/*for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++)
			for(uint k=0;k<3;k++){
				ex z=pow(fmasses[1][i][i]/McH,2); mixes.M[tQuark][iUp][i][i]/MR,2);
				llgammaH0M[j][k]=llgammaH0M[j][k]+(mixes.VN[0][1][j][i].conjugate()*mixes.VN[0][1][k][i]+mixes.VN[0][1][i][j]*mixes.VN[0][1][i][k].conjugate())/pow(mixes.M[tLepton][iDown][i][i],2)*(2*z+6*z*z*log(z))/6;
				llgammaH0M[j][k]=llgammaH0M[j][k]+(mixes.VN[0][1][i][j]*mixes.VN[0][1][k][i])/mixes.M[tLepton][iDown][i][i]/mixes.M[tLepton][iDown][j][j]*(3*z+2*z*log(z));
				
				llgammaH0E[j][k]=llgammaH0E[j][k]+(mixes.VN[0][1][j][i].conjugate()*mixes.VN[0][1][k][i]-mixes.VN[0][1][i][j]*mixes.VN[0][1][i][k].conjugate())/pow(mixes.M[tLepton][iDown][i][i],2)*(2*z+6*z*z*log(z))/6;
				llgammaH0E[j][k]=llgammaH0E[j][k]+(mixes.VN[0][1][i][j]*mixes.VN[0][1][k][i])/mixes.M[tLepton][iDown][i][i]/mixes.M[tLepton][iDown][j][j]*(3*z+2*z*log(z));
				}
				*/
	llgammaH0M=Matrix(g*g/(16*Pi*4*Pi*MW*MW*McH*McH))*mixes.M[tLepton][iDown]*llgammaH0M;
	llgammaH0E=Matrix(g*g/(16*Pi*4*Pi*MW*MW*McH*McH))*mixes.M[tLepton][iDown]*llgammaH0E;
	
	Matrix llgamma, llgamma2;
	
	for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++){
		//if(j<i) llgamma[i][j]=(llgammaCH[i][j]*llgammaCH[i][j].conjugate()+llgammaH0E[i][j]*llgammaH0E[i][j].conjugate()+llgammaH0M[i][j]*llgammaH0M[i][j].conjugate())*g*g*(1-cos2)*pow((pow(mixes.M[tLepton][iDown][i][i],2)-pow(mixes.M[tLepton][iDown][j][j],2))/mixes.M[tLepton][iDown][i][i],3)/(4*Pi);
	    ex mmuon=mixes.M[tLepton][iDown][i][i];
			ex A,B;
			
	    if(j<i){ for(uint k=0;k<3;k++){
			  ex mtau=mixes.M[tLepton][iDown][k][k];
			  B+=-mixes.VN[tLepton][1][k][j].conjugate()*mixes.VN[tLepton][1][k][i]/(12*pow(McH,2));
			  B+=mixes.N[tLepton][1][k][j].conjugate()*mixes.N[tLepton][1][k][i]/12*(pow(MR,-2)+pow(MI,-2));
			  
			  A+=mixes.N[tLepton][1][j][k]*mixes.N[tLepton][1][i][k].conjugate()*(pow(MR,-2)+pow(MI,-2))/12;
			  A+=mixes.N[tLepton][1][j][k]*mixes.N[tLepton][1][k][i]/mtau/mmuon*(Fh2(pow(mtau/MR,2))-Fh2(pow(mtau/MI,2)))/4;
			}
			 llgamma[i][j]=(A*A.conjugate()+B*B.conjugate())*alpha*pow(mmuon,5)*GF*GF/(128*pow(Pi,4));
			 }
		else if(j==i){
			for(uint k=0;k<3;k++){
			  ex mtau=mixes.M[tLepton][iDown][k][k];
			  B+=-mixes.VN[tLepton][1][k][j].conjugate()*mixes.VN[tLepton][1][k][i]/(12*pow(McH,2));
			  B+=mixes.N[tLepton][1][k][j].conjugate()*mixes.N[tLepton][1][k][i]/12*(pow(MR,-2)+pow(MI,-2));
			  B+=mixes.N[tLepton][1][j][k].conjugate()*mixes.N[tLepton][1][i][k]/12*(pow(MR,-2)+pow(MI,-2));
			  }
		    llgamma[i][j]=-B*GF*sqrt(1/2)/(8*pow(Pi,2))*2*mmuon; //e (GeV)^-1=1/(51e6)(e cm) where e=sqrt(alpha*4*Pi)
		  }
		}
	add("mutoegamma",llgamma[1][0],new limitedobs(planck/2.197034e-6*2.4e-12));
	add("tautoegamma",llgamma[2][0],new limitedobs(planck/290.6e-15*3.3e-8));
	add("tautomugamma",llgamma[2][1],new limitedobs(planck/290.6e-15*4.4e-8));
	
	/*add("d_e",abs(llgamma[0][0].imag_part()),new limitedobs(10.5e-28*51e6));
	add("d_mu",abs(llgamma[1][1].imag_part()),new limitedobs(1.9e-19*51e6));
	add("d_tau",llgamma[2][2].imag_part(),new gaussobs(-0.85e-17*51e6,0.825/0.85));
	cout<<"EDM: "<<llgamma[0][0].subs(conjtoabs).subs(replacements).imag_part()<<endl;
	add("a_mu",-llgamma[1][1].real_part()*2*mixes.M[tLepton][iDown][1][1],new gaussobs(3e-9,1.0/3.0));
	*/
	/*llgammaCH=Matrix(g*g/(16*Pi*4*Pi*MW*MW*McH*McH*12))*mixes.M[tQuark][iDown]*mixes.VN[1][1].conjugate()*mixes.VN[1][1]; //4+1
	//Matrix llgammaH0M,llgammaH0E;
	for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++)
			for(uint k=0;k<3;k++){
				ex z=pow(mixes.M[tQuark][iUp][i][i]/MR,2);
				llgammaH0M[j][k]=llgammaH0M[j][k]+(mixes.VN[1][1][j][i].conjugate()*mixes.VN[1][1][k][i]+mixes.VN[1][1][i][j]*mixes.VN[1][1][i][k].conjugate())/pow(mixes.M[tQuark][iDown][i][i],2)*(2*z+6*z*z*log(z))/6;
				llgammaH0M[j][k]=llgammaH0M[j][k]+(mixes.VN[1][1][i][j]*mixes.VN[1][1][k][i])/mixes.M[tQuark][iDown][i][i]/mixes.M[tQuark][iDown][j][j]*(3*z+2*z*log(z));
				
				llgammaH0E[j][k]=llgammaH0E[j][k]+(mixes.VN[1][1][j][i].conjugate()*mixes.VN[1][1][k][i]-mixes.VN[1][1][i][j]*mixes.VN[1][1][i][k].conjugate())/pow(mixes.M[tQuark][iDown][i][i],2)*(2*z+6*z*z*log(z))/6;
				llgammaH0E[j][k]=llgammaH0E[j][k]+(mixes.VN[1][1][i][j]*mixes.VN[1][1][k][i])/mixes.M[tQuark][iDown][i][i]/mixes.M[tQuark][iDown][j][j]*(3*z+2*z*log(z));
				}
	llgammaH0M=Matrix(g*g/(16*Pi*4*Pi*MW*MW*McH*McH))*mixes.M[tQuark][iDown]*llgammaH0M;
	llgammaH0E=Matrix(g*g/(16*Pi*4*Pi*MW*MW*McH*McH))*mixes.M[tQuark][iDown]*llgammaH0E;
	
	//Matrix llgamma;
	for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++) 
		if(j<i) {llgamma[i][j]=(llgammaCH[i][j]*llgammaCH[i][j].conjugate()+llgammaH0E[i][j]*llgammaH0E[i][j].conjugate()+llgammaH0M[i][j]*llgammaH0M[i][j].conjugate())*g*g*(1-cos2)*pow((pow(mixes.M[tQuark][iDown][i][i],2)-pow(mixes.M[tQuark][iDown][j][j],2))/mixes.M[tQuark][iDown][i][i],3)/(4*Pi);
				//llgamma[i][j]=llgamma[i][j].subs(lst(abs(wild()*pow(MH0,-2))==abs(wild())*pow(MH0,-2)));
				}
	    else llgamma[i][j]=0;
	 */

	
	push_back(prediction(new calcubtosgamma2(mixes)));
	    
	//add("btosgamma",llgamma[2][1],new gaussobs(3.55e-4,sqrt(2)*0.25/3.55),1);
	    //cout<<csrc<<llgamma[2][1]<<endl;
	    //cout<<latex;

	
	BR_Htotaunu=(CHdecaycoupling(chiggs,tau,neutrino)+3*CHdecaycoupling(chiggs,strange,charm))/factor(CHdecaycoupling(chiggs,Fermion(tLepton,iDown),neutrino)+3*CHdecaycoupling(chiggs,Fermion(tQuark,iDown),charm)+3*CHdecaycoupling(chiggs,Fermion(tQuark,iDown),up));
	BR_Htotaunu=BR_Htotaunu.subs(replacements);
	
	//BR_toptoHq=decaywidth(top,bottom,chiggs);
	//ex toptoWb=decaywidth(top,bottom,wboson);
	//BR_toptoHq=BR_toptoHq/(BR_toptoHq+toptoWb);
	//BR_toptoHq=BR_toptoHq.subs(replacements);
	
	//cout<<"toptoWb "<<toptoWb.subs(replacements).evalf()<<endl;
	
	//b to c tau- nu/b to c e- nu
	//ex btocR=decaywidth(bottom,charm,tau,neutrino,sVector)/(decaywidth(bottom,charm,electron,neutrino,sVector)+decaywidth(bottom,charm,muon,neutrino,sVector));
	//cout<<btocR.subs(replacements)<<endl;
	
	ex BtoDtaunu,BtoD2taunu, BtoDtaunuSM, KtoPi;
	for(uint i=0; i<3; i++){
		ex Wcoup=wboson.couplingL(charm,bottom)*wboson.couplingdaggerL(tau,Fermion(tLepton,iUp,FFlavour(i)));
		if(Wcoup.subs(replacements)==ex(0)) continue;
		ex chcoup_Wcoup=-pow(MW/McH,2)*(chiggs.couplingR(charm,bottom)+chiggs.couplingL(charm,bottom))*chiggs.couplingdaggerL(tau,Fermion(tLepton,iUp,FFlavour(i)))/Wcoup;
		ex chcoup2_Wcoup=-pow(MW/McH,2)*(chiggs.couplingR(charm,bottom)-chiggs.couplingL(charm,bottom))*chiggs.couplingdaggerL(tau,Fermion(tLepton,iUp,FFlavour(i)))/Wcoup;
		
		BtoDtaunuSM+=Wcoup*Wcoup.conjugate();
		BtoDtaunu+=Wcoup*Wcoup.conjugate()*(1+1.5*chcoup_Wcoup.real_part()+chcoup_Wcoup.conjugate()*chcoup_Wcoup);
		BtoD2taunu+=Wcoup*Wcoup.conjugate()*(1+0.12*chcoup2_Wcoup.real_part()+0.05*chcoup2_Wcoup.conjugate()*chcoup2_Wcoup);	
	}
	lst r2(pow(mixes.V[1][1][2].imag_part(),2)==pow(abs(mixes.V[1][1][2]),2)-pow(mixes.V[1][1][2].real_part(),2));
	r2.append(pow(mixes.V[0][2][2].imag_part(),2)==pow(abs(mixes.V[0][2][2]),2)-pow(mixes.V[0][2][2].real_part(),2));
	
	r2.append(mixes.M[1][0][1][1]==0);
	r2.append(pow(abs(mixes.V[0][2][2]),2)==1-pow(abs(mixes.V[0][1][2]),2)-pow(abs(mixes.V[0][0][2]),2));
	r2.append(pow(abs(mixes.V[0][2][1]),2)==1-pow(abs(mixes.V[0][1][1]),2)-pow(abs(mixes.V[0][0][1]),2));
	r2.append(abs(sqrt(ex(2))* GF)==sqrt(ex(2))* GF);
	
	BtoDtaunuSM=collect_common_factors(BtoDtaunuSM.subs(conjtoabs).subs(r2));
	BtoDtaunu=collect_common_factors(BtoDtaunu.subs(conjtoabs).subs(r2));
	
	BtoDtaunuR=(BtoDtaunu/BtoDtaunuSM).subs(replacements).real_part();
	
	BtoD2taunu=BtoD2taunu.subs(conjtoabs).subs(r2);
	BtoD2taunuR=(BtoD2taunu/BtoDtaunuSM).subs(replacements).real_part();
	
	
	//cout<<"BtoDtaunu/BtoDtaunuSM "<<expand(BtoDtaunu/BtoDtaunuSM)<<endl;
	iBDtaunu=size();
	add("BtoDtaunu_BtoDtaunuSM",BtoDtaunu/BtoDtaunuSM,new gaussobs(440.0/296, 1.4*58.0/440));
	
	iBD2taunu=size();
	//cout<<"BtoD2taunu/BtoD2taunuSM "<<1+collect_common_factors(expand(BtoD2taunu/BtoDtaunuSM-1))<<endl;
	add("BtoD2taunu_BtoD2taunuSM",BtoD2taunu/BtoDtaunuSM,new gaussobs(332.0/252, 1.4*24.0/332.0));
	
	
	
	for(uint j=0; j<2; j++){
		ex KtoPimunu, KtoPimunuSM;
	for(uint i=0; i<3; i++){
		ex Wcoup=wboson.couplingL(up,strange)*wboson.couplingdaggerL(Fermion(tLepton,iDown,FFlavour(j)),Fermion(tLepton,iUp,FFlavour(i)));
		if(Wcoup.subs(replacements)==ex(0)) continue;
		ex chcoup_Wcoup=-pow(MW/McH,2)*(chiggs.couplingR(up,strange)+chiggs.couplingL(up,strange))\
				*chiggs.couplingdaggerL(muon,Fermion(tLepton,iUp,FFlavour(i)))/Wcoup*(pow(MKp,2)-pow(Mpip,2))\
				/(mixes.mass(Fermion(tLepton,iDown,FFlavour(j)))*(mixes.mass(strange)-mixes.mass(up)));
		chcoup_Wcoup=collect_common_factors(expand(chcoup_Wcoup));
		KtoPimunuSM+=collect_common_factors(expand(Wcoup*Wcoup.conjugate()));
		KtoPimunu+=collect_common_factors(expand(Wcoup*Wcoup.conjugate()*pow(1+chcoup_Wcoup,2)));
	}
	    KtoPimunuSM=collect_common_factors(expand(KtoPimunuSM.subs(conjtoabs).subs(r2)));
	    KtoPimunu=collect_common_factors(expand(KtoPimunu.subs(conjtoabs).subs(r2)));
	    KtoPimunu=expand(KtoPimunu.subs(replacements).real_part().subs(lst(abs(wild()*pow(MR,-2))==abs(wild())*pow(MR,-2))).subs(lst(log(wild()*pow(MR,-2))==log(wild())-2*log(MR))));
	    KtoPimunu=expand(KtoPimunu.evalf());
	    KtoPimunuSM=expand(KtoPimunuSM.subs(replacements).real_part().subs(lst(abs(wild()*pow(MR,-2))==abs(wild())*pow(MR,-2))).subs(lst(log(wild()*pow(MR,-2))==log(wild())-2*log(MR))));
	    KtoPimunuSM=expand(KtoPimunuSM.evalf());
		KtoPi+=0.5*log(KtoPimunu/KtoPimunuSM);
	}
	
	add("KtoPi",KtoPi/(pow(MKp,2)-pow(Mpip,2)),new gaussobs(0.08, 0.11/0.08));
	
	
	//add("b to c tau- nu/b to c e- nu", decaywidth(bottom,charm,electron,neutrino), new limitedobs(planck/290.6e-15*2.7e-8));

    double fD=0.207;
    ex DDbar=ex(std::pow(fD,2))*mesonmixing(MD0,charm,up);
    DDbar=expand(DDbar.subs(replacements).subs(lst(abs(wild()*pow(MR,-2))==abs(wild())*pow(MR,-2))).subs(lst(log(wild()*pow(MR,-2))==log(wild())-2*log(MR))));
	DDbar=expand(DDbar.evalf());
    ex aDDbar=sqrt(DDbar.real_part()*DDbar.real_part()+DDbar.imag_part()*DDbar.imag_part());   
    add("DDbar",aDDbar,new limitedobs(9.47e-15));
    cout<<DDbar<<endl;
//2|M12|<6.6e-15GeV     
	
    double fK=0.156;
    ex KKbar=ex(std::pow(fK,2))*mesonmixing(MK0,strange,down);
    KKbar=expand(KKbar.subs(replacements).subs(lst(abs(wild()*pow(MR,-2))==abs(wild())*pow(MR,-2))).subs(lst(log(wild()*pow(MR,-2))==log(wild())-2*log(MR))));
	KKbar=expand(KKbar.evalf());
	ex aKKbar=sqrt(KKbar.real_part()*KKbar.real_part()+KKbar.imag_part()*KKbar.imag_part());   
    add("KKbar",aKKbar,new limitedobs(3.5e-15));
    ex eK=0.94*imag_part(KKbar)/3.5e-15/sqrt(2);
    //add("a_eK",abs(eK),new limitedobs(2.2e-3));
    add("a_eK",abs(eK),new limitedobs(20*0.0114e-3));
    cout<<abs(KKbar)<<endl;
    
    double fB=0.189;
    ex Vtb=mixes.V[tQuark][2][2]/mixes.V[tQuark][2][2].conjugate();
    ex Vtd=mixes.V[tQuark][2][0]/mixes.V[tQuark][2][0].conjugate();
    ex Vts=mixes.V[tQuark][2][1]/mixes.V[tQuark][2][1].conjugate();
    
    ex BBbar=1+ex(std::pow(fB,2))*mesonmixing(MB0,bottom,down)/(3.337e-13*Vtb*Vtd.conjugate());
    add("BBbarimag",imag_part(BBbar),new gauss2obs(-0.199,0.062));
	add("BBbarreal",real_part(BBbar),new gauss2obs(0.823,0.143));
    cout<<BBbar<<endl;
    BBbar=3.337e-13*Vtb*Vtd.conjugate();
    cout<<"Bbar "<<(abs(imag_part(BBbar))/abs(BBbar)).subs(replacements)<<endl;
	double fBs=0.225;
    ex BsBsbar=1+ex(std::pow(fBs,2))*mesonmixing(MBs0,bottom,strange)/(1.186e-11*Vtb*Vts.conjugate());
    add("BsBsbarimag",imag_part(BsBsbar),new gauss2obs(0,0.1));
	add("BsBsbarreal",real_part(BsBsbar),new gauss2obs(0.965,0.133));
    cout<<BsBsbar<<endl;
    BsBsbar=1.186e-11*Vtb*Vts.conjugate();
    cout<<"Bbar "<<(abs(imag_part(BsBsbar))/abs(BsBsbar)).subs(replacements)<<endl;
	
    ex McH2=McH*McH;
    ex MR2=MR*MR;
    ex MI2=MI*MI;
    
    ex cu=collect_common_factors(expand(chiggs.couplingL(top,bottom)))/mixes.mass(top)/(g/MW/sqrt(ex(2)))/mixes.V[1][2][2];
    cout<<"cu "<<cu<<endl;
    ex Zbb=(cu-0.72)/McH;
    add("Zbb",Zbb,new limitedobs(0.0024));
    cout<<"Zbb "<<Zbb<<endl;
    cout<<"SIZE "<<size()<<endl;
      
   push_back(prediction(new calcuOblique()));
}

~BGL(){
	delete cBmumu;
	delete cBsmumu;
	}

ex Y(ex x) const{
		return 1.0113*x/8/(1-x)*(4-x+3*x*log(x)/(1-x));
		}
	
ex GW(ex x) const{
	return (94-x*(179+x*(-55+12*x)))/(36*pow(1-x,3))+x*(16+x*(-32+9*x))*log(x)/(6*pow(1-x,4));
	}

ex GH1(ex x) const{
	return x*(x*((39-14*x)*x-6)+6*x*(3*x-8)*log(x)-19)/(36*pow(x-1,4));
	//return -x/12;
	}

ex GH2(ex x) const{
	return x*((x-1)*(11*x-21)+(16-6*x)*log(x))/(6*pow(x-1,3));
	//return x/2;
	}

ex FW(ex x) const{
	return (94-x*(179+x*(-55+12*x)))/(36*pow(1-x,3))+x*(16+x*(-32+9*x))*log(x)/(6*pow(1-x,4));
	}

ex FH1(ex x) const{
	return -x/12;
	}

ex FH2(ex x) const{
	return x/2;
	}
	
ex Fh1(ex x) const{
	//return (2*x+3*pow(x,2)-6*pow(x,3)+pow(x,4)+6*pow(x,2)*log(x))/(6*pow(1-x,4));
	return x/3;
	}

ex Fh2(ex x) const{
	//return (-3*x+4*pow(x,2)-pow(x,3)-2*x*log(x))/pow(1-x,3);
	return -2*(3/2+log(x))*x;
	}

ex A0(ex x) const{
	return x*(2+3*x-6*x*x+ x*x*x+6*x*log(x))/(24*pow(1-x,4));
	}

ex A1(ex x) const{
	return x*(-3+4*x-x*x-2*log(x))/(4*pow(1-x,3));
	}

ex A2(ex x) const{
	return x/(6*pow(1-x,3))*((-7+5*x-8*x*x)/6.0+x*log(x)/(1-x)*(-2+3*x));
	}

ex A3(ex x) const{
	return (-3+8*x-5*x*x+(6*x-4)*log(x))*x/(6*pow(1-x,3));
	}
		
void add(const char * s, ex pred, observable * ob, bool sb=0){
	//cout<<s<<endl;
	//cout<<"prediction symb"<<pred<<endl;
	//,pow(sin(wild()), 2) == 1-pow(cos(wild()), 2)
	//ex p=expand(pred.subs(replacements).real_part().subs(lst(abs(wild()*pow(MR,-2))==abs(wild())*pow(MR,-2))).subs(lst(log(wild()*pow(MR,-2))==log(wild())-2*log(MR))));
	
	ex p=pred.subs(replacements).real_part();
	p=collect_common_factors(expand(p.evalf()));
	FUNCP_CUBA fp;
	
	lst l(tanb,McH,MR,MI);	
	
	for(uint i=0;i<3;i++){
		l.append(Mu[i]);
		l.append(Md[i]);
	}
	if(sb) push_back(prediction(ob,p));
	else {
	compile_ex(lst(p), l, fp);
	//cout<<"prediction numeric"<<p<<endl;
	//cout<<"exp "<<ob->expected()<<endl<<endl;
	push_back(prediction(ob,fp));
   }
	}

int veto(const parameters & p, int max=0) const{
	if(!p.isvalid()) return 1;
	if(max==1){
	double mr=p[1].value+p[2].value;
	if(mr<10 || mr>10000) return 1;
	mr+=p[3].value;
	if(mr<10 || mr>10000) return 1;
	return 0;
	}
	else{
	double mr=p[1].value+p[2].value;
	if(mr<10 || mr>mmmax) return 1;
	mr+=p[3].value;
	if(mr<10 || mr>mmmax) return 1;
	return 0;	
	}
	}

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
	
	return pp;
}

double bsgammawidth(double tanb,double McH,double MR,double MI, int option=0){
	parameters p=generateparameters();
	p[0].value=pow(10.0,tanb);
	p[1].value=McH;
	p[2].value=MR;
	p[3].value=MI;
	
	calcubtosgamma2 cal(mixes);
	
	return cal.width(p,option);
}
/*
ex decaywidth2(const Fermion& f1, const Fermion& ff2, const Fermion& ff3, const Fermion& ff4, BSpin s=sAny) const{
	
	Fermion f2=ff2,f3=ff3, f4=ff4;
	
	ex ret=0;
	
	realsymbol q1("q1"), q2("q2"), q3("q3"), q4("q2"), s3("s3");
	ex s2=pow(mixes.mass(f1),2);
	
	for(uint j=fElectron;j<=fTau;j++) 
	if(ff2.flavour==fAny || ff2.flavour==j){
		f2.flavour=(FFlavour)j;
	for(uint k=fElectron;k<=fTau;k++) 
	if(ff3.flavour==fAny || ff3.flavour==k){
		f3.flavour=(FFlavour)k;
	for(uint l=fElectron;l<=fTau;l++)
	if(ff4.flavour==fAny || ff4.flavour==l){
		f4.flavour=(FFlavour)l;		
		ex v1=0, v2=0;
		ex mq1=mixes.mass(f1),mq2=mixes.mass(f2),mq3=mixes.mass(f3),mq4=mixes.mass(f4);
		ex m2q1=mq1*mq1, m2q2=mq2*mq2, m2q3=mq3*mq3, m2q4=mq4*mq4;
		ex s4=m2q1+m2q2+m2q3+m2q4-s2-s3;
		
		scalar_products sp;
		sp.add(q2, q3, (s4-m2q2-m2q3)/2); 
		sp.add(q4, q3, (s2-m2q4-m2q3)/2); 
		sp.add(q2, q4, (s3-m2q2-m2q4)/2); 
  
		sp.add(q2, q2, m2q2);
		sp.add(q3, q3, m2q3);
		sp.add(q4, q4, m2q4);

		ex q10=(s2+mq1*mq1-mq2*mq2)/(2*sqrt(s2)), lq1l=sqrt(q10*q10-mq1*mq1);
		ex q30=(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
    
	for(uint i=0;i<bosons.size();i++)if(bosons[i].s==s || s==sAny){
		if(bosons[i].s==0){
			ex a=-(bosons[i].couplingR(f2,f1)-bosons[i].couplingL(f2,f1))*bosons[i].couplingdaggerL(f3,f4)*s2/(mq1+mq2)/pow(bosons[i].mass,2);
			v1=v1+a*dirac_gammaL();
			v2=v2+a.conjugate()*dirac_gammaR();
			a=-(bosons[i].couplingR(f2,f1)-bosons[i].couplingL(f2,f1))*bosons[i].couplingdaggerR(f3,f4)*s2/(mq1+mq2)/pow(bosons[i].mass,2);
			v1=v1+a*dirac_gammaR();
			v2=v2+a.conjugate()*dirac_gammaL();
		}
		else{
			ex sl=(dirac_slash(q3,4)+dirac_slash(q4,4));
			ex a=(bosons[i].couplingR(f2,f1)-bosons[i].couplingL(f2,f1))*bosons[i].couplingdaggerL(f3,f4)/pow(bosons[i].mass,2);
			v1=v1+a*sl*dirac_gammaL();
			v2=v2+a.conjugate()*sl*dirac_gammaL();
			a=(bosons[i].couplingR(f2,f1)-bosons[i].couplingL(f2,f1))*bosons[i].couplingdaggerR(f3,f4)/pow(bosons[i].mass,2);
			v1=v1+a*sl*dirac_gammaR();
			v2=v2+a.conjugate()*sl*dirac_gammaR();
		}
	}
	ex vq3=dirac_slash(q3,4)-mq3*dirac_ONE(), vq4=dirac_slash(q4,4)+mq4*dirac_ONE();
	ex dt=dirac_trace(vq3*v1*vq4*v2).simplify_indexed(sp);
	cout<<"dt: "<<dt<<endl;
	ex result=expand(dt*4*lq3l/s2/Pi/128);

	ret+=result;
	}
	}
	
	return collect_common_factors(ret.subs(conjtoabs));
	//return expand(ret.subs(lst(exp(-I*wild())==1/exp(I*wild()),sin(wild())==sqrt(1-pow(cos(wild()),2)))));
}
*/

ex decaywidth(const Fermion& ff1, const Fermion& ff2, const Fermion& ff3, const Fermion& ff4, BSpin s=sAny) const{
	multivector<ex,4> a(0,bosons.size(),2,2,2);
	vector<ex> mass(bosons.size(),0);
	vector<int> op(bosons.size(),0);
	ex ret=0;
	Fermion f1=ff1, f2=ff2,f3=ff3, f4=ff4;
	
	
	for(uint i=fElectron;i<=fTau;i=i+1) 
	if(ff1.flavour==fAny || ff1.flavour==i){
		f1.flavour=(FFlavour)i;
	for(uint j=fElectron;j<=fTau;j++)
	if(ff2.flavour==fAny || ff2.flavour==j){
		f2.flavour=(FFlavour)j;
	for(uint k=fElectron;k<=fTau;k++) 
	if(ff3.flavour==fAny || ff3.flavour==k){
		f3.flavour=(FFlavour)k;
	for(uint l=fElectron;l<=fTau;l++)
	if(ff4.flavour==fAny || ff4.flavour==l){
		f4.flavour=(FFlavour)l;		
	for(uint i=0;i<bosons.size();i++) if(bosons[i].s==s || s==sAny){
			op[i]=bosons[i].s;
			mass[i]=bosons[i].mass;
			a[i][0][0][0]=bosons[i].couplingdaggerL(f2,f1)*bosons[i].couplingL(f3,f4);
			a[i][0][0][1]=bosons[i].couplingdaggerL(f2,f1)*bosons[i].couplingR(f3,f4);
			a[i][0][1][0]=bosons[i].couplingdaggerR(f2,f1)*bosons[i].couplingL(f3,f4);
			a[i][0][1][1]=bosons[i].couplingdaggerR(f2,f1)*bosons[i].couplingR(f3,f4);
			
			a[i][1][0][0]=bosons[i].couplingdaggerL(f3,f1)*bosons[i].couplingL(f2,f4);
			a[i][1][0][1]=bosons[i].couplingdaggerL(f3,f1)*bosons[i].couplingR(f2,f4);
			a[i][1][1][0]=bosons[i].couplingdaggerR(f3,f1)*bosons[i].couplingL(f2,f4);
			a[i][1][1][1]=bosons[i].couplingdaggerR(f3,f1)*bosons[i].couplingR(f2,f4);		
	}
	
	ret+=wc.get_integral_symb(a,mass,op,mixes.mass(f1));
	//ret+=wc.get_integral(a,mass,op,mixes.massnum(f1),mixes.massnum(f2),mixes.massnum(f3),mixes.massnum(f4))/pow(mixes.massnum(f1),5)*pow(mixes.mass(f1),5);	
	}}}}
	if(ff2.flavour==ff4.flavour) ret=ret/2;
	return collect_common_factors(ret.subs(conjtoabs));
	//return expand(ret.subs(lst(exp(-I*wild())==1/exp(I*wild()),sin(wild())==sqrt(1-pow(cos(wild()),2)))));
}

ex get_integral_symb(const multivector<ex,3>& a, ex m1) const{
	realsymbol s2("s2"), s3("s3");
	realsymbol q1("q1"), q2("q2"), q3("q3"), q4("q4");
			
	ex m2q1=m1*m1;
	
	ex vq1=dirac_slash(q2,4)+dirac_slash(q3,4)+dirac_slash(q4,4)+m1*dirac_ONE(), vq2=dirac_slash(q2,4);
	ex vq3=dirac_slash(q3,4), vq4=dirac_slash(q4,4);
	
	ex s4=m2q1-s2-s3;
    scalar_products sp;
    sp.add(q2, q3, (s4)/2); 
    sp.add(q4, q3, (s2)/2); 
    sp.add(q2, q4, (s3)/2); 
  
    sp.add(q2, q2, 0);
    sp.add(q3, q3, 0);
    sp.add(q4, q4, 0);
    
    multivector<ex,2> v(0,2,2);
	v[0][0]=dirac_gammaL(); v[0][1]=dirac_gammaR();
	v[1][0]=dirac_gammaR(); v[1][1]=dirac_gammaL();
	
	multivector<ex,5> traces(0,2,2,2,2,2);
	    for(uint k=0;k<2;k++)
		  for(uint l=0;l<2;l++)
		  for(uint m=0;m<2;m++)
				for(uint n=0;n<2;n++){
				ex vk=v[k][0];
				ex vm=v[m][0];
				ex vl=v[l][1];
				ex vn=v[n][1];
	
				traces[k][l][m][n][0]=dirac_trace(vq2*vk*vq1*vl)*dirac_trace(vq3*vm*vq4*vn);
				traces[k][l][m][n][1]=-dirac_trace(vq2*vk*vq1*vl*vq3*vm*vq4*vn);
			}
			
	    for(uint k=0;k<2;k++)
			for(uint l=0;l<2;l++)
				for(uint m=0;m<2;m++)
				for(uint n=0;n<2;n++)
						for(uint o=0;o<2;o++)
					{
							traces[k][l][m][n][o]=(traces[k][l][m][n][o]).simplify_indexed(sp);
						}
							
	ex q10=(s2+m1*m1)/(2*sqrt(s2)), lq1l=(m1*m1-s2)/(2*sqrt(s2));
    ex q30=sqrt(s2)/2, lq3l=q30;
    ex q20=(m1*m1-s2)/(2*m1), lq2l=q20;
    
    ex total=0;
    for(uint k=0;k<2;k++)
	for(uint l=0;l<2;l++)
	for(uint m=0;m<2;m++)
	for(uint n=0;n<2;n++)
	for(uint r=0;r<2;r++)
	for(uint s=0;s<2;s++){
		ex coup=a[r][k][m]*a[s][l][n].conjugate();
		ex integrand=traces[k][l][m][n][(r+s)%2];
		integrand=expand(integral(s3, 0, m1*m1-s2, integrand).eval_integ()/lq1l/sqrt(s2)*lq2l/m1/m1);
		//double mm2=0, mm3=0, m4=0;
		ex result=integral(s2,0,m1*m1,integrand).eval_integ()/pow(Pi,3)/512;
	    ex partial=result*coup;
		total=total+partial;
		}
	return total;
}

ex decaywidthtest2(const Fermion& ff1) const{
	multivector<ex,3> a(0,2,2,2);
	symbol gLL("g_{LL}"),gLR("g_{LR}"),gRL("g_{RL}"),gRR("g_{RR}"),cLL("c_{LL}"),cLR("c_{LR}"),cRL("c_{RL}"),cRR("c_{RR}");
		
			a[0][0][0]=gLL;
			a[0][0][1]=gLR;
			a[0][1][0]=gRL;
			a[0][1][1]=gRR;
			
			a[1][0][0]=cLL;
			a[1][0][1]=cLR;
			a[1][1][0]=cRL;
			a[1][1][1]=cRR;		
	
	ex ret=get_integral_symb(a,mixes.mass(ff1));
	//ret+=wc.get_integral(a,mass,op,mixes.massnum(f1),mixes.massnum(f2),mixes.massnum(f3),mixes.massnum(f4))/pow(mixes.massnum(f1),5)*pow(mixes.mass(f1),5);	
	
	return collect_common_factors(ret.subs(conjtoabs));
	//return expand(ret.subs(lst(exp(-I*wild())==1/exp(I*wild()),sin(wild())==sqrt(1-pow(cos(wild()),2)))));
}

/*
ex decaywidthtest(const Fermion& f1, const Fermion& f2, const Fermion& ff3, const Fermion& ff4, BSpin s=sAny) const{
	
	Fermion f1=ff1, f2=ff2,f3=ff3, f4=ff4;
	
	ex ret=0;
	
	realsymbol q1("q1"), q2("q2"), q3("q3"), q4("q4");
	symbol gL("gL"), gR("gR");
	ex s2("s2"), s3("s3");
	
	for(uint i=fElectron;i<=fTau;i=i+1) 
	if(ff1.flavour==fAny || ff1.flavour==i){
		f1.flavour=(FFlavour)i;
	for(uint j=fElectron;j<=fTau;j++) 
	if(ff2.flavour==fAny || ff2.flavour==j){
		f2.flavour=(FFlavour)j;	
	for(uint k=fElectron;k<=fTau;k++) 
	if(ff3.flavour==fAny || ff3.flavour==k){
		f3.flavour=(FFlavour)k;
	for(uint l=fElectron;l<=fTau;l++)
	if(ff4.flavour==fAny || ff4.flavour==l){
		f4.flavour=(FFlavour)l;
		ex v1=0, v2=0;
		ex mq1=mixes.mass(f1),mq2=mixes.mass(f2),mq3=mixes.mass(f3),mq4=mixes.mass(f4);
		ex m2q1=mq1*mq1, m2q2=mq2*mq2, m2q3=mq3*mq3, m2q4=mq4*mq4;
		
		ex s4=m2q1+m2q2+m2q3+m2q4-s2-s3;
	
		scalar_products sp;
			
		sp.add(q4, q3, (s2-m2q4-m2q3)/2);
		sp.add(q3, q3, m2q3);
		sp.add(q4, q4, m2q4);
		ex q10=(s2+mq1*mq1-mq2*mq2)/(2*sqrt(s2)), lq1l=sqrt(q10*q10-mq1*mq1);
		ex q30=(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
		
		ex vq1=dirac_slash(q2,4)+dirac_slash(q3,4)+dirac_slash(q4,4)+mq1*dirac_ONE(), vq2=dirac_slash(q2,4)+mq2*dirac_ONE();
 	    ex vq3=dirac_slash(q3,4)+mq3*dirac_ONE(), vq4=dirac_slash(q4,4)-mq4*dirac_ONE();
		
		ex a;
		
		a=gL;
		v1=v1+a*dirac_gammaL();
		v2=v2+a.conjugate()*dirac_gammaR();
		a=gR;
		v1=v1+a*dirac_gammaR();
		v2=v2+a.conjugate()*dirac_gammaL();
		
		/ *}
		else{
			ex sl=(dirac_slash(q3,4)+dirac_slash(q4,4));
			ex a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingL(f3,f4)/pow(bosons[i].mass,2);
			v1=v1+a*sl*dirac_gammaL();
			v2=v2+a.conjugate()*sl*dirac_gammaL();
			a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingR(f3,f4)/pow(bosons[i].mass,2);
			v1=v1+a*sl*dirac_gammaR();
			v2=v2+a.conjugate()*sl*dirac_gammaR();
		}* /
		
	ex vq3=dirac_slash(q3,4)+mq3*dirac_ONE(), vq4=dirac_slash(q4,4)-mq4*dirac_ONE();
	ex dt=dirac_trace(vq3*v1*vq4*v2).simplify_indexed(sp);
	ex result=expand(dt*4*lq3l/s2/Pi/128);

	ret+=result;
	}
	}
	lst ltest;
	ltest.append(conjugate(gL)==pow(abs(gL),2)/gL);
	ltest.append(conjugate(gR)==pow(abs(gR),2)/gR);
	return pow(meson.decay_factor,2)*collect_common_factors(ret.subs(conjtoabs).subs(ltest));
	//return expand(ret.subs(lst(exp(-I*wild())==1/exp(I*wild()),sin(wild())==sqrt(1-pow(cos(wild()),2)))));
}
*/

ex gRR2(const Fermion& f1, const Fermion& f3) const{
	
	ex ret1=0,ret2=0;
	Fermion f2(tLepton,iUp);
	Fermion f4(tLepton,iUp);
	
	for(uint k=fElectron;k<=fTau;k++){ 
		f2.flavour=(FFlavour)k;
	for(uint l=fElectron;l<=fTau;l++){
		f4.flavour=(FFlavour)l;		
	for(uint i=0;i<bosons.size();i++)
		if(bosons[i].s==sScalar) {
			ex x=bosons[i].couplingdaggerR(f2,f1)*bosons[i].couplingL(f3,f4)/pow(bosons[i].mass,2);
			ret1+=x*x.conjugate();
		}
		else if(bosons[i].s==sVector) {
			ex x=bosons[i].couplingdaggerL(f2,f1)*bosons[i].couplingL(f3,f4)/pow(bosons[i].mass,2);
			ret2+=x*x.conjugate();
		}
	}}
	//r2.append();
	ret2=ret2.subs(conjtoabs);
	ret1=ret1.subs(conjtoabs);
	for(uint i=0;i<3;i++){
		ret1=collect_common_factors(expand(ret1.subs(pow(abs(mixes.V[0][2][i]),2)==1-pow(abs(mixes.V[0][1][i]),2)-pow(abs(mixes.V[0][0][i]),2))));
		ret2=collect_common_factors(expand(ret2.subs(pow(abs(mixes.V[0][2][i]),2)==1-pow(abs(mixes.V[0][1][i]),2)-pow(abs(mixes.V[0][0][i]),2))));
	}

	cout<<ret2<<endl;
	return collect_common_factors(ret1/ret2);
}

ex tautomu_tautoe() const{
	
	ex ret1=0,ret2=0, rety1=0, rety2=0;
	
	Fermion f1(tLepton,iDown,fTau);
	Fermion f31(tLepton,iDown,fMuon);
	Fermion f32(tLepton,iDown,fElectron);
	
	Fermion f2(tLepton,iUp);
	Fermion f4(tLepton,iUp);
	
	
	for(uint k=fElectron;k<=fTau;k++){ 
		f2.flavour=(FFlavour)k;
	for(uint l=fElectron;l<=fTau;l++){
		f4.flavour=(FFlavour)l;
		ex x1=0, x2=0, y1=0, y2=0;		
	    for(uint i=0;i<bosons.size();i++){
		if(bosons[i].s==sScalar) {
			x1+=bosons[i].couplingdaggerR(f2,f1)*bosons[i].couplingL(f31,f4)/pow(bosons[i].mass,2);
			x2+=bosons[i].couplingdaggerR(f2,f1)*bosons[i].couplingL(f32,f4)/pow(bosons[i].mass,2);
		}
		else if(bosons[i].s==sVector) {
			y1+=bosons[i].couplingdaggerL(f2,f1)*bosons[i].couplingL(f31,f4)/pow(bosons[i].mass,2);
			y2+=bosons[i].couplingdaggerL(f2,f1)*bosons[i].couplingL(f32,f4)/pow(bosons[i].mass,2);
		}
		}
		ret1+=(x1*y1.conjugate()).real_part();
		ret2+=(x2*y2.conjugate()).real_part();
		rety1+=y1*y1.conjugate();
		rety2+=y2*y2.conjugate();
	}}
	ret2=(ret2/rety2*mixes.mass(f32)/mixes.mass(f1)).subs(conjtoabs);
	ret1=(ret1/rety1*mixes.mass(f31)/mixes.mass(f1)).subs(conjtoabs);
	for(uint i=0;i<3;i++){
		ret1=collect_common_factors(expand(ret1.subs(pow(abs(mixes.V[0][2][i]),2)==1-pow(abs(mixes.V[0][1][i]),2)-pow(abs(mixes.V[0][0][i]),2))));
		ret2=collect_common_factors(expand(ret2.subs(pow(abs(mixes.V[0][2][i]),2)==1-pow(abs(mixes.V[0][1][i]),2)-pow(abs(mixes.V[0][0][i]),2))));
	}
	
	ex x=pow(mixes.mass(f31)/mixes.mass(f1),2);
	ex F1=1-8*x+8*pow(x,3)-pow(x,4)-12*pow(x,2)*log(x);
	ex g1=1+9*x-9*pow(x,2)-pow(x,3)+6*x*(1+x)*log(x);
	ex N1=1+gRR2(f1,f31)/4;
	
	x=pow(mixes.mass(f32)/mixes.mass(f1),2);
	
	ex F2=1-8*x+8*pow(x,3)-pow(x,4)-12*pow(x,2)*log(x);
	ex g2=1+9*x-9*pow(x,2)-pow(x,3)+6*x*(1+x)*log(x);
	ex N2=1+gRR2(f1,f32)/4;
	
	return collect_common_factors(N1*(F1+2/N1*ret1*g1)/N2/(F2+2/N2*ret2*g2)*F2/F1);
}

ex mesondw(const Meson & meson, const Fermion& ff3, const Fermion& ff4, BSpin s=sAny) const{
	
	const Fermion& f1(meson.q1), f2(meson.q2);
	ex mesonmass=meson.mass;
	
	Fermion f3=ff3, f4=ff4;
	
	ex ret=0;
	
	realsymbol q3("q3"), q4("q4");
	ex s2=pow(mesonmass,2);
		
	for(uint k=fElectron;k<=fTau;k++) 
	if(ff3.flavour==fAny || ff3.flavour==k){
		f3.flavour=(FFlavour)k;
	for(uint l=fElectron;l<=fTau;l++)
	if(ff4.flavour==fAny || ff4.flavour==l){
		f4.flavour=(FFlavour)l;
		ex v1=0, v2=0;
		ex mq1=mixes.mass(f1),mq2=mixes.mass(f2),mq3=mixes.mass(f3),mq4=mixes.mass(f4);
		ex m2q1=mq1*mq1, m2q2=mq2*mq2, m2q3=mq3*mq3, m2q4=mq4*mq4;
		scalar_products sp;
		sp.add(q4, q3, (s2-m2q4-m2q3)/2);
		sp.add(q3, q3, m2q3);
		sp.add(q4, q4, m2q4);
		ex q10=(s2+mq1*mq1-mq2*mq2)/(2*sqrt(s2)), lq1l=sqrt(q10*q10-mq1*mq1);
		ex q30=(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
    
	for(uint i=0;i<bosons.size();i++) if(bosons[i].s==s || s==sAny){
		if(bosons[i].s==0){
			ex a=-(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingL(f3,f4)*s2/(mq1+mq2)/pow(bosons[i].mass,2);
			v1=v1+a*dirac_gammaL();
			v2=v2+a.conjugate()*dirac_gammaR();
			a=-(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingR(f3,f4)*s2/(mq1+mq2)/pow(bosons[i].mass,2);
			v1=v1+a*dirac_gammaR();
			v2=v2+a.conjugate()*dirac_gammaL();
		}
		else{
			ex sl=(dirac_slash(q3,4)+dirac_slash(q4,4));
			ex a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingL(f3,f4)/pow(bosons[i].mass,2);
			v1=v1+a*sl*dirac_gammaL();
			v2=v2+a.conjugate()*sl*dirac_gammaL();
			a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingR(f3,f4)/pow(bosons[i].mass,2);
			v1=v1+a*sl*dirac_gammaR();
			v2=v2+a.conjugate()*sl*dirac_gammaR();
		}
	}
	ex vq3=dirac_slash(q3,4)+mq3*dirac_ONE(), vq4=dirac_slash(q4,4)-mq4*dirac_ONE();
	ex dt=dirac_trace(vq3*v1*vq4*v2).simplify_indexed(sp);
	ex result=expand(dt*4*lq3l/s2/Pi/128);

	ret+=result;
	}
	}
	
	return pow(meson.decay_factor,2)*collect_common_factors(ret.subs(conjtoabs));
	//return expand(ret.subs(lst(exp(-I*wild())==1/exp(I*wild()),sin(wild())==sqrt(1-pow(cos(wild()),2)))));
}


ex mesondwtest(const Meson & meson, const Fermion& ff3, const Fermion& ff4, BSpin s=sAny) const{
	
	const Fermion& f1(meson.q1), f2(meson.q2);
	ex mesonmass=meson.mass;
	
	Fermion f3=ff3, f4=ff4;
	
	ex ret=0;
	
	realsymbol q3("q3"), q4("q4");
	symbol gL("gL"), gR("gR"),gVL("gVL"), gVR("gVR");
	symbol gS("gS"), gP("gP"), gA("gA");
	
	ex s2=pow(mesonmass,2);
		
	for(uint k=fElectron;k<=fTau;k++) 
	if(ff3.flavour==fAny || ff3.flavour==k){
		f3.flavour=(FFlavour)k;
	for(uint l=fElectron;l<=fTau;l++)
	if(ff4.flavour==fAny || ff4.flavour==l){
		f4.flavour=(FFlavour)l;
		ex v1=0, v2=0;
		ex mq1=mixes.mass(f1),mq2=mixes.mass(f2),mq3=mixes.mass(f3),mq4=mixes.mass(f4);
		ex m2q1=mq1*mq1, m2q2=mq2*mq2, m2q3=mq3*mq3, m2q4=mq4*mq4;
		scalar_products sp;
		sp.add(q4, q3, (s2-m2q4-m2q3)/2);
		sp.add(q3, q3, m2q3);
		sp.add(q4, q4, m2q4);
		ex q10=(s2+mq1*mq1-mq2*mq2)/(2*sqrt(s2)), lq1l=sqrt(q10*q10-mq1*mq1);
		ex q30=(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
		
		ex a;
	/*	a=-gL*s2/(mq1+mq2);
		v1=v1+a*dirac_gammaL();
		v2=v2+a.conjugate()*dirac_gammaR();
		a=-gR*s2/(mq1+mq2);
		v1=v1+a*dirac_gammaR();
		v2=v2+a.conjugate()*dirac_gammaL();
		
		ex sl=(dirac_slash(q3,4)+dirac_slash(q4,4));
		a=gA;
		v1=v1+a*sl*dirac_gamma5();
		v2=v2+a.conjugate()*sl*dirac_gamma5();
*/
		a=-gS*s2/(mq1+mq2);
		v1=v1+a*dirac_ONE();	
		v2=v2+a.conjugate()*dirac_ONE();
		a=-gP*s2/(mq1+mq2);
		v1=v1+a*dirac_gamma5();
		v2=v2-a.conjugate()*dirac_gamma5();
		ex sl=(dirac_slash(q3,4)+dirac_slash(q4,4));
	//	a=gA;
	//	v1=v1+a*sl*dirac_gamma5();
	//	v2=v2+a.conjugate()*sl*dirac_gamma5();
		
		/*}
		else{
			ex sl=(dirac_slash(q3,4)+dirac_slash(q4,4));
			ex a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingL(f3,f4)/pow(bosons[i].mass,2);
			v1=v1+a*sl*dirac_gammaL();
			v2=v2+a.conjugate()*sl*dirac_gammaL();
			a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingR(f3,f4)/pow(bosons[i].mass,2);
			v1=v1+a*sl*dirac_gammaR();
			v2=v2+a.conjugate()*sl*dirac_gammaR();
		}*/
		
	ex vq3=dirac_slash(q3,4)+mq3*dirac_ONE(), vq4=dirac_slash(q4,4)-mq4*dirac_ONE();
	ex dt=dirac_trace(vq3*v1*vq4*v2).simplify_indexed(sp);
	ex result=expand(dt*4*lq3l/s2/Pi/128);

	ret+=result;
	}
	}
	lst ltest;
	ltest.append(conjugate(gL)==pow(abs(gL),2)/gL);
	ltest.append(conjugate(gR)==pow(abs(gR),2)/gR);
	ltest.append(conjugate(gS)==pow(abs(gS),2)/gS);
	ltest.append(conjugate(gP)==pow(abs(gP),2)/gP);
	ltest.append(conjugate(gA)==pow(abs(gA),2)/gA);
	
	return pow(meson.decay_factor,2)*collect_common_factors(expand(ret.subs(conjtoabs).subs(ltest)));
	//return expand(ret.subs(lst(exp(-I*wild())==1/exp(I*wild()),sin(wild())==sqrt(1-pow(cos(wild()),2)))));
}

ex fermiontomeson(const Fermion& ff4, const Fermion& ff3, const Meson & meson, BSpin s=sAny) const{
	
	const Fermion& f1(meson.q1), f2(meson.q2);
	ex mesonmass=meson.mass;
	
	Fermion f3=ff3, f4=ff4;
	
	ex ret=0;
	
	realsymbol q3("q3"), q4("q4");
	ex s2=pow(mesonmass,2);
		
	for(uint k=fElectron;k<=fTau;k++) 
	if(ff3.flavour==fAny || ff3.flavour==k){
		f3.flavour=(FFlavour)k;
	for(uint l=fElectron;l<=fTau;l++)
	if(ff4.flavour==fAny || ff4.flavour==l){
		f4.flavour=(FFlavour)l;
		ex v1=0, v2=0;
		ex mq1=mixes.mass(f1),mq2=mixes.mass(f2),mq3=mixes.mass(f3),mq4=mixes.mass(f4);
		ex m2q1=mq1*mq1, m2q2=mq2*mq2, m2q3=mq3*mq3, m2q4=mq4*mq4;
		scalar_products sp;
		sp.add(q4, q3, -(s2-m2q4-m2q3)/2);
		sp.add(q3, q3, m2q3);
		sp.add(q4, q4, m2q4);
		//ex qm0=(s2-mq3*mq3+mq4*mq4)/(2*mq4), lqml=sqrt(qm0*qm0-s2);
		ex q30=(-s2+mq3*mq3+mq4*mq4)/(2*mq4), lq3l=sqrt(q30*q30-mq3*mq3);
        //ex q30=-(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
    
	for(uint i=0;i<bosons.size();i++)if(bosons[i].s==s || s==sAny){
		if(bosons[i].s==0){
			ex a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingL(f3,f4)*s2/(mq1+mq2)/pow(bosons[i].mass,2);
			v1=v1+a*dirac_gammaL();
			v2=v2+a.conjugate()*dirac_gammaR();
			a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingR(f3,f4)*s2/(mq1+mq2)/pow(bosons[i].mass,2);
			v1=v1+a*dirac_gammaR();
			v2=v2+a.conjugate()*dirac_gammaL();
		}
		else{
			ex sl=(dirac_slash(q3,4)+dirac_slash(q4,4));
			ex a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingL(f3,f4)/pow(bosons[i].mass,2);
			v1=v1+a*sl*dirac_gammaL();
			v2=v2+a.conjugate()*sl*dirac_gammaL();
			a=(bosons[i].couplingdaggerR(f2,f1)-bosons[i].couplingdaggerL(f2,f1))*bosons[i].couplingR(f3,f4)/pow(bosons[i].mass,2);
			v1=v1+a*sl*dirac_gammaR();
			v2=v2+a.conjugate()*sl*dirac_gammaR();
		}
	}
	ex vq3=dirac_slash(q3,4)+mq3*dirac_ONE(), vq4=dirac_slash(q4,4)+mq4*dirac_ONE();
	ex dt=dirac_trace(vq3*v1*vq4*v2).simplify_indexed(sp);
	ex result=expand(dt*2*lq3l/mq4/mq4/Pi/128);

	ret+=result;
	}
	}
	
	
	
	return pow(meson.decay_factor,2)*collect_common_factors(ret.subs(conjtoabs));
	//return expand(ret.subs(lst(exp(-I*wild())==1/exp(I*wild()),sin(wild())==sqrt(1-pow(cos(wild()),2)))));
}

ex fermiontomesontest(const Fermion& ff4, const Fermion& ff3, const Meson & meson, BSpin s=sAny) const{
	
	const Fermion& f1(meson.q1), f2(meson.q2);
	ex mesonmass=meson.mass;
	
	Fermion f3=ff3, f4=ff4;
	
	ex ret=0;
	
	realsymbol q3("q3"), q4("q4");
	
	symbol sL("sL"), sR("sR"), vL("vL"), vR("vR");
	ex s2=pow(mesonmass,2);
		
	for(uint k=fElectron;k<=fTau;k++) 
	if(ff3.flavour==fAny || ff3.flavour==k){
		f3.flavour=(FFlavour)k;
	for(uint l=fElectron;l<=fTau;l++)
	if(ff4.flavour==fAny || ff4.flavour==l){
		f4.flavour=(FFlavour)l;
		ex v1=0, v2=0;
		ex mq1=mixes.mass(f1),mq2=mixes.mass(f2),mq3=mixes.mass(f3),mq4=mixes.mass(f4);
		ex m2q1=mq1*mq1, m2q2=mq2*mq2, m2q3=mq3*mq3, m2q4=mq4*mq4;
		scalar_products sp;
		sp.add(q4, q3, -(s2-m2q4-m2q3)/2);
		sp.add(q3, q3, m2q3);
		sp.add(q4, q4, m2q4);
		//ex qm0=(s2-mq3*mq3+mq4*mq4)/(2*mq4), lqml=sqrt(qm0*qm0-s2);
		ex q30=(-s2+mq3*mq3+mq4*mq4)/(2*mq4), lq3l=sqrt(q30*q30-mq3*mq3);
        //ex q30=-(s2+mq3*mq3-mq4*mq4)/(2*sqrt(s2)), lq3l=sqrt(q30*q30-mq3*mq3);
    
	
			ex a=sL;
			v1=v1+a*dirac_gammaL();
			v2=v2+a.conjugate()*dirac_gammaR();
			a=sR;
			v1=v1+a*dirac_gammaR();
			v2=v2+a.conjugate()*dirac_gammaL();

			ex sl=(dirac_slash(q3,4)+dirac_slash(q4,4));
			a=vL;
			v1=v1+a*sl*dirac_gammaL();
			v2=v2+a.conjugate()*sl*dirac_gammaL();
			a=vR;
			v1=v1+a*sl*dirac_gammaR();
			v2=v2+a.conjugate()*sl*dirac_gammaR();
			
	ex vq3=dirac_slash(q3,4)+mq3*dirac_ONE(), vq4=dirac_slash(q4,4)+mq4*dirac_ONE();
	ex dt=dirac_trace(vq3*v1*vq4*v2).simplify_indexed(sp);
	ex result=expand(dt*2*lq3l/mq4/mq4/Pi/128);

	ret+=result;
	}
	}
	
	return pow(meson.decay_factor,2)*collect_common_factors(ret.subs(conjtoabs));
	//return expand(ret.subs(lst(exp(-I*wild())==1/exp(I*wild()),sin(wild())==sqrt(1-pow(cos(wild()),2)))));
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

ex CHdecaycoupling(Boson higgs, const Fermion& ff3, const Fermion& ff4) const{
	
	Fermion f3=ff3, f4=ff4;
	ex ret=0;
	for(uint k=fElectron;k<=fTau;k++) 
	if(ff3.flavour==fAny || ff3.flavour==k){
		f3.flavour=(FFlavour)k;
	for(uint l=fElectron;l<=fTau;l++)
	if(ff4.flavour==fAny || ff4.flavour==l){
		f4.flavour=(FFlavour)l;		
		ret+=higgs.couplingdaggerL(f3,f4)*higgs.couplingdaggerL(f3,f4).conjugate()+higgs.couplingdaggerR(f3,f4)*higgs.couplingdaggerR(f3,f4).conjugate();
	}}
	return collect_common_factors(ret.subs(conjtoabs));
}


double BranchingRatio(double * xx,double * p){
	return ex_to<numeric>(BR_Htotaunu.subs(tanb==pow(10.0,xx[0])).evalf()).to_double();
}


double topBranchingRatio(double * xx,double * p){
	return ex_to<numeric>(BR_toptoHq.subs(lst(tanb==pow(10.0,xx[0]),McH==xx[1])).evalf()).to_double();
}

widthcalc wc;

const double planck;
const possymbol GF, MZ, MW, Mh;
const constant Mpip, Mpi0, MBp,MB0,MBs0, MKp,MK0,MDp,MD0,MDsp,MDs0;
const constant Fpi, FB,FBs, FK,FD,FDs;
ex cos2, g, alpha;
const possymbol tanb, cp, McH, MR, MI, rho;
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

double mmmax, stepsize;

calcuBmumu * cBmumu;
calcuBmumu * cBsmumu;


};


}
#endif
