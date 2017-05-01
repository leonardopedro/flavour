#include <iostream>
#include <fstream>

using namespace std;
typedef unsigned int uint;
int main(){

string g[3]={"0","1","2"};
string u[3]={"0","1"};

string sg[3]={"1st","2nd","3rd"};
string su[3]={"Down","Up"};

for(uint qup=0;qup<2;qup++)
for(uint gL=0;gL<3;gL++)
{
for(uint lup=0;lup<2;lup++)
for(uint gQ=0;gQ<3;gQ++){
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
  double McH=0,MR=0,MI=0;
  ff>>McH>>MR>>MI;
  cout<<name<<" "<<McH<<" "<<MR<<" "<<MI<<endl;
  ff.close();
}
cout<<endl;
}
return 0;
}
