#include <iostream>
#include <fstream>


using namespace std;
typedef unsigned int uint;
int main(){

string g[3]={"0","1","2"};
string u[3]={"0","1"};

string sg[3]={"1st","2nd","3rd"};
string su[3]={"Down","Up"};

for(uint qup=0;qup<2;qup++){
cout<<"\\pagebreak\n";
cout<<"\\begin{figure}[!htb]\n\\centering"<<endl;
for(uint gL=0;gL<3;gL++)
for(uint lup=0;lup<2;lup++){
cout<<"\\begin{tikzpicture}[every node/.style={anchor=south west,inner sep=0pt},"<<endl;
cout<<"x=0.307692307692\\textwidth,y=0.230769230769\\textwidth,]"<<endl;	

for(uint gQ=0;gQ<3;gQ++)
{
		char name[5]="0000";
		name[0]+=gL;
		name[1]+=gQ;
		name[2]+=lup;
		name[3]+=qup;
  
  ifstream ff((string("pdfs/T_BD/maxs_")+string(name)+string(".out")).c_str());
  if(!ff.is_open()){
		cout<<"ERROR: maxs_*.out not found"<<endl;
		return 1;
	}
  int tx=0, ty=0;
  double xi=0,xxi=0.26,xw=0.384615384615,yi=0.25;
  if(gQ>0) {xi=0.24+gQ; xxi+=gQ; tx=192; xw=0.307692307692;}
  if(gL<2 || lup<1) {ty=144; yi=0.01;}
  
  double McH=0,MR=0,MI=0;
  ff>>McH>>MR>>MI;
  ff.close();
  int mcH(McH+0.5),mR(MR+0.5),mI(MI+0.5);

  cout<<"\\node at ("<<xi<<",0) {\\includegraphics[trim="<<tx<<" "<<ty<<" 0 0,clip,"<<endl;
  cout<<"width="<<xw<<"\\textwidth]{../pdfs/T_BD/pdf_"<<name<<".png}};"<<endl;
  cout<<"\\node at ("<<xxi<<","<<yi<<") {\\begin{minipage}{0.1\\textwidth}{"<<endl;
  cout<<"\\scriptsize \\begin{align*} M_{H^+}&>"<<mcH<<"\\\\[-4pt]"<<endl;
  cout<<"M_{R^0}&>"<<mR<<"\\\\[-4pt]"<<endl;
  cout<<"M_{I^0}&>"<<mI<<" (GeV)\\end{align*}}\\end{minipage}};"<<endl;
  
  
}
cout<<"\\end{tikzpicture}\\\\"<<endl;
}
cout<<"\\caption{BGL Quarks: gen ";
cout<<"; Leptons: gen FCNC chargedlepton/neutrino on the left/right.}"<<endl;
cout<<"\\end{figure}"<<endl<<endl;
}
 
 
 for(uint qup=0;qup<2;qup++)
 for(uint lup=0;lup<2;lup++)
 for(uint gL=0;gL<3;gL++)
 for(uint gQ=0;gQ<3;gQ++){
  cout<<"export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH && ./teste "<<gL<<" "<<gQ<<" "<<lup<<" "<<qup;
  cout<<" > "<<"teste"<<gL<<gQ<<lup<<qup<<".out"<<endl;
 }
 
 for(uint lup=0;lup<2;lup++)
 for(uint qup=0;qup<2;qup++)
 for(uint gL=0;gL<3;gL++)
 for(uint gQ=0;gQ<3;gQ++){
  cout<<"./update "<<gL<<" "<<gQ<<" "<<lup<<" "<<qup;
  cout<<" > "<<"teste"<<gL<<gQ<<lup<<qup<<".out"<<endl;
 }
return 0;
}
