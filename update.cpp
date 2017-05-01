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
#include <iostream>
#include <fstream>

using namespace std;

/**
* @brief the main function takes the arguments file1 file2 output, merges the resuts in ROOT file1 and file2 producing the output ROOT file 
*/
int main(int argc, char* argv[]){
    // Check the number of parameters
    
    if(argc < 4){
        std::cerr<<"Usage: "<<argv[0]<<" file1 file2 output"<<std::endl;
        return 1;}
    
	TH2F *limits4,*limits_tb_MR;
	TH2F *limits_MR_MI,*limits_MR_McH,*limits_MI_McH;
	TH2F *limits4_,*limits_tb_MR_;
	TH2F *limits_MR_MI_,*limits_MR_McH_,*limits_MI_McH_;
	TVectorD *vllmax=NULL,*vllmax_=NULL;
	
	uint npoints=200;
	
	TFile *f=new TFile(argv[2],"update");
	if(!f->IsOpen()) cout<<"NOFILE"<<endl;
	f->GetObject("vllmax;1",vllmax_);
	f->GetObject("limits4;1",limits4_);
	f->GetObject("limits_tb_MR;1",limits_tb_MR_);
	f->GetObject("limits_MR_MI;1",limits_MR_MI_);
	f->GetObject("limits_MR_McH;1",limits_MR_McH_);
	f->GetObject("limits_MI_McH;1",limits_MI_McH_);
	
	if(!vllmax_) cout<<"ERROR"<<endl;
	
	TFile *f2=new TFile(argv[1],"update");
	if(!f2->IsOpen()) cout<<"NOFILE"<<endl;

	f2->GetObject("vllmax;1",vllmax);
	f2->GetObject("limits4;1",limits4);
	f2->GetObject("limits_tb_MR;1",limits_tb_MR);
	f2->GetObject("limits_MR_MI;1",limits_MR_MI);
	f2->GetObject("limits_MR_McH;1",limits_MR_McH);
	f2->GetObject("limits_MI_McH;1",limits_MI_McH);
	
	double c=(*vllmax_)[0];
	if(c>(*vllmax)[0]) (*vllmax)[0]=c;
	for(uint i=0;i<npoints;i++)
		for(uint j=0;j<npoints;j++) {
			c=limits4_->GetBinContent(i+1,j+1);
			if(c>limits4->GetBinContent(i+1,j+1)) limits4->SetBinContent(i+1,j+1,c);
			c=limits_tb_MR_->GetBinContent(i+1,j+1);
			if(c>limits_tb_MR->GetBinContent(i+1,j+1)) limits_tb_MR->SetBinContent(i+1,j+1,c);
			c=limits_MR_MI_->GetBinContent(i+1,j+1);
			if(c>limits_MR_MI->GetBinContent(i+1,j+1)) limits_MR_MI->SetBinContent(i+1,j+1,c);
			c=limits_MR_McH_->GetBinContent(i+1,j+1);
			if(c>limits_MR_McH->GetBinContent(i+1,j+1)) limits_MR_McH->SetBinContent(i+1,j+1,c);
			c=limits_MI_McH_->GetBinContent(i+1,j+1);
			if(c>limits_MI_McH->GetBinContent(i+1,j+1)) limits_MI_McH->SetBinContent(i+1,j+1,c);
			}
	TFile *f3=new TFile(argv[3],"recreate");
	if(!f3->IsOpen()) cout<<"NOFILE"<<endl;
	
	vllmax->Write("vllmax");		
	limits4->Write();
	limits_MR_MI->Write();
    limits_MR_McH->Write();
    limits_MI_McH->Write();
    limits_tb_MR->Write();
    
    f3->Close();
	f2->Close();
	f->Close();
	

	return 0;
	
}

