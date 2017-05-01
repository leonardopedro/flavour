#include <iostream>

using namespace std;
typedef unsigned int uint;
int main(){

string g[3]={"0","1","2"};
string u[3]={"0","1"};

string sg[3]={"1st","2nd","3rd"};
string su[3]={"Down","Up"};

for(uint qup=0;qup<2;qup++)
for(uint gQ=0;gQ<3;gQ++)
{
cout<<"\\pagebreak\n";
for(uint gL=0;gL<3;gL++){
cout<<"\\begin{figure}[!htb]\n\\centering"<<endl;
for(uint lup=0;lup<2;lup++){
cout<<"\\includegraphics[width=0.49\\textwidth]{../pdfs/T_BD/pdf_";
    cout<<g[gL]<<g[gQ]<<u[lup]<<u[qup];
    cout<<".png}"<<endl;
}
cout<<"\\caption{BGL Quarks: gen "<<sg[gQ]<<" FCNC "<<su[qup];
cout<<"; Leptons: gen "<<sg[gL]<<" FCNC chargedlepton/neutrino on the left/right.}"<<endl;
cout<<"\\end{figure}"<<endl<<endl;
}
}
return 0;
}
