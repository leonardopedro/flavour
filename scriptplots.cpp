#include <string>
#include <iostream>
#include <fstream>

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

using namespace std;

int main(){

string g[3]={"0","1","2"};
string u[3]={"0","1"};

string sg[3]={"1st","2nd","3rd"};
string su[3]={"Down","Up"};

ofstream f[2];

f[0].open("draft/large_fig1.tex");
f[1].open("draft/large_fig2.tex");

for(uint qup=0;qup<2;qup++){
	
for(uint gL=0;gL<3;gL++)
for(uint lup=0;lup<2;lup++){
f[qup]<<"\\begin{tikzpicture}[every node/.style={anchor=south west,inner sep=0pt},"<<endl;
f[qup]<<"x=0.307692307692\\textwidth,y=0.230769230769\\textwidth,]"<<endl;	

for(uint gQ=0;gQ<3;gQ++)
{
		char name[5]="0000";
		name[0]+=gL;
		name[1]+=gQ;
		name[2]+=lup;
		name[3]+=qup;
		char name_[8]="0 0 0 0";
		name_[0]+=gL;
		name_[2]+=gQ;
		name_[4]+=lup;
		name_[6]+=qup;

  //system((string("./update pdfs/T_BD4/h")+string(name)+string(".root pdfs/T_BD3/h")+string(name)+string(".root pdfs/T_BD/h")+string(name)+string(".root")).c_str());
  
  system((string("cd pdfs/T_BD; ../../draw ./h")+string(name)+string(".root ")+string(name_)).c_str());
  ifstream ff((string("pdfs/T_BD/maxs_")+string(name)+string(".out")).c_str());
  if(!ff.is_open()){
		cout<<"ERROR: maxs_*.out not found"<<endl;
		return 1;
	}
  int tx=0, ty=0;
  double xi=0,xxi=0.26,xw=0.384615384615,yi=0.25;
  yi=0.25;
  if(gQ>0) {xi=0.24+gQ; xxi+=gQ; tx=192; xw=0.307692307692;}
  if(gL<2 || lup<1) {ty=144; yi=0.01;}
  
  double McH=0,MR=0,MI=0,llmax=-1000,tbmax=0,McHmax=1000,MRmax=1000,MImax=1000;
  double eK__=0,sm__=0,charged__=0,neutral__=0,neutralR__=0,neutralI__=0,all__=0;
  ff>>McH>>MR>>MI;
  if(gQ>0 && qup) {
	  ff>>McH>>MR>>MI;
	  ff>>McH>>MR>>MI;
  }else{ ff>>llmax>>llmax>>llmax;
    ff>>llmax>>llmax>>llmax;
  }
  ff>>llmax>>tbmax>>McHmax>>MRmax>>MImax;
  //ff>>eK__;
  //ff>>sm__>>charged__>>neutral__>>neutralR__>>neutralI__>>all__;
  ff.close();
  int eK_(eK__*100+0.5);
  int sm_(sm__*100+0.5);
  int charged_(charged__*100+0.5);
  int neutral_(neutral__*100+0.5);
  int neutralR_(neutralR__*100+0.5);
  int neutralI_(neutralI__*100+0.5);
  int all_(all__*100+0.5);
  int mcH(McH+0.5),mR(MR+0.5),mI(MI+0.5);
  double tmax(int(tbmax*100+0.5)/100.0);
  int mcHmax(McHmax+0.5),mRmax(MRmax+0.5),mImax(MImax+0.5);
  f[qup]<<"\\node at ("<<xi<<",0) {\\includegraphics[trim="<<tx<<" "<<ty<<" 0 0,clip,"<<endl;
  f[qup]<<"width="<<xw<<"\\textwidth]{../pdfs/T_BD/pdf_"<<name<<".png}};"<<endl;
  f[qup]<<"\\node at ("<<xxi<<","<<yi<<") {\\begin{minipage}{0.1\\textwidth}{"<<endl;
  f[qup]<<"\\scriptsize \\begin{align*}"<<endl;
  //f[qup]<<tmax<<" & \\ "<<mcHmax<<"\\ "<<mRmax<<"\\ "<<mImax<<"\\\\[-4pt]"<<endl;
  //f[qup]<<sm_<<" & \\ "<<charged_<<"\\ "<<neutral_<<"\\ "<<neutralR_<<"\\ "<<neutralI_<<"\\ "<<all_<<"\\\\[-4pt]"<<endl;
  //f[qup]<<" & \\ "<<eK_<<"\\\\[-4pt]"<<endl;
  f[qup]<<mcH<<"&<M_{H^+}/\\GeV\\\\[-4pt]"<<endl;
  f[qup]<<mR<<"&<M_{R^0}/\\GeV\\\\[-4pt]"<<endl;
  f[qup]<<mI<<"&<M_{I^0}/\\GeV\\end{align*}}\\end{minipage}};"<<endl;
  
  
}
f[qup]<<"\\end{tikzpicture}\\\\"<<endl;
}
f[qup].close();
}
 
return 0;
}
