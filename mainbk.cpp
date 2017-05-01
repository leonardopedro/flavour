#include "MCMC.h"
#include "BGL.h"
#include "TF2.h"
#include "TProfile3D.h"
#include "THStack.h"
#include <cln/cln.h>
#include <cln/float.h>
int main(){
	Digits=5;
	cln::cl_inhibit_floating_point_underflow=1;
	
	
	
	
	//TH1F * pdf1=new TH1F("pdf1","pdf1",npoints,10,500);
	
	
	//TGraph * chi2=new TGraph(npoints);
	double lmax=0;
	double gaussmax=0;
	double tanbmax=0,mcHmax=0, Btaunu=0,BDtaunu=0,BD2taunu=0;
	uint gLm=0,gQm=0,lupm=0,qupm=0;
	
	for(uint gL=0;gL<3;gL++)
	for(uint gQ=0;gQ<3;gQ++)
	for(uint lup=0;lup<2;lup++)
	for(uint qup=0;qup<2;qup++){
		BGL* m=new BGL(gL,gQ,lup,qup);
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
	uint npoints=100;	
	//TH2F * pdf=new TH2F("pdf","pdf",npoints,-7/log(10.0),7/log(10.0),npoints,10,500);
	double init=-3, final=3;
	double init2=1, final2=4;
	
	uint steps=1e4;
	Proposal prop1(m,m->iBtaunu);
	prop1.findPeaks();
	Proposal prop2(m,m->iBDtaunu);
	prop2.findPeaks();
	Proposal prop3(m,m->iBD2taunu);
	prop3.findPeaks();
	Proposal prop4(m,-1);
	prop4.findPeaks();
	
	tanbmax=pow(10.0,prop4.floatPeak.pr[0].value);
	mcHmax=pow(10.0,prop4.floatPeak.pr[1].value);
	
	/*
	
	TH2F * limits1=new TH2F("limits1","Likelihood",npoints,init,final,npoints,init2,final2);
	TH2F * limits2=new TH2F("limits2","Likelihood",npoints,init,final,npoints,init2,final2);
	TH2F * limits3=new TH2F("limits3","Likelihood",npoints,init,final,npoints,init2,final2);
	TH2F * limits4=new TH2F("limits4","Likelihood",npoints,init,final,npoints,init2,final2);
	TH2F * limits5=new TH2F("limits5","Likelihood",npoints,init,final,npoints,init2,final2);
	TProfile2D * like1=new TProfile2D("like1","like",npoints,init,final,npoints,init2,final2);
	TProfile2D * like2=new TProfile2D("like2","like",npoints,init,final,npoints,init2,final2);
	TProfile2D * like3=new TProfile2D("like3","like",npoints,init,final,npoints,init2,final2);
	TProfile2D * like4=new TProfile2D("like4","like",npoints,init,final,npoints,init2,final2);
	THStack hs("hs","test stacked histograms");
	
	
	for(uint i=steps;i;i--){
		prop1.getNextPoint();
		prop2.getNextPoint();
		prop3.getNextPoint();
		prop4.getNextPoint();
		
		double Btaunu=m->loglike(prop1.floatPeak.pr,m->iBtaunu);
		double BDtaunu=m->loglike(prop2.floatPeak.pr,m->iBDtaunu);
		double BD2taunu=m->loglike(prop3.floatPeak.pr,m->iBD2taunu);
		double total=m->loglike(prop4.floatPeak.pr,-1);
		
		like1->Fill(prop1.floatPeak.pr[0].value, prop1.floatPeak.pr[1].value,Btaunu);
		like2->Fill(prop2.floatPeak.pr[0].value, prop2.floatPeak.pr[1].value,BDtaunu);
		like3->Fill(prop3.floatPeak.pr[0].value, prop3.floatPeak.pr[1].value,BD2taunu);
		like4->Fill(prop4.floatPeak.pr[0].value, prop4.floatPeak.pr[1].value,total);
					
		if(i%(steps/10)==0) cout<<"Steps "<<i<<endl;
		//pdf1->Fill(prop.floatPeak.pr[1].value);
	}
		
	for(uint i=0;i<npoints;i++)
		for(uint j=0;j<npoints;j++)
		{
			double bcmax=0;
			int binmax=like1->GetBin(i+1,j+1);
			//for(uint k=0;k<npoints;k++){
			//int binmax=like1->GetBin(i+1,j+1,j+1);
			//if(like1->GetBinEntries(bin)){
			//	double bc=like1->GetBinContent(bin);
			//	if(bc>bcmax || binmax==-1) {bcmax=bc; binmax=bin;}
			//	}
			//}
			double Btaunu=0;
			double BDtaunu=0;
			double BD2taunu=0;
			double rest=0;
			double sum=0;
			
			if(like1->GetBinEntries(binmax)) {
				Btaunu=like1->GetBinContent(binmax);
				//cout<<i<<" "<<j<<" "<<Btaunu<<" ";
				sum+=Btaunu;
				Btaunu=TMath::Prob(-2*Btaunu,1);
				//cout<<Btaunu<<endl;
				if(Btaunu<0.05) Btaunu=0;
			}
			if(like2->GetBinEntries(binmax)) {
				BDtaunu=like2->GetBinContent(binmax);
				sum+=BDtaunu;
				BDtaunu=TMath::Prob(-2*BDtaunu,1);
				if(BDtaunu<0.05) BDtaunu=0;
			}
			if(like3->GetBinEntries(binmax)) {
				BD2taunu=like3->GetBinContent(binmax);
				sum+=BD2taunu;
				BD2taunu=TMath::Prob(-2*BD2taunu,1);
				if(BD2taunu<0.05) BD2taunu=0;
			}
			if((like3->GetBinEntries(binmax) && like2->GetBinEntries(binmax) && like1->GetBinEntries(binmax) && sum>-20)){
				sum=TMath::Prob(-2*sum,3);
				if(sum>lmax){
				lmax=sum;
				tanbmax=pow(10.0,(i+0.5)*1.0/npoints*6-3);
				mcHmax=pow(10.0,(j+0.5)*1.0/npoints*3+1);
			}
			}else sum=0;
			
			if(like4->GetBinEntries(binmax)){
				rest=like4->GetBinContent(binmax);
				rest=TMath::Prob(-2*rest,m->size());
			}
			int scale=100;
			limits1->SetBinContent(i+1,j+1,Btaunu*scale);
			limits2->SetBinContent(i+1,j+1,BDtaunu*scale);
			limits3->SetBinContent(i+1,j+1,BD2taunu*scale);
			limits4->SetBinContent(i+1,j+1,rest*scale);
			limits5->SetBinContent(i+1,j+1,sum*scale);
		}
	
	limits3->GetXaxis()->SetTitle("tan\\beta*MW/MH");
    limits3->GetYaxis()->SetTitle("cotan\\beta*MW/MH");
    limits1->SetStats(0);
    limits2->SetStats(0);
    limits3->SetStats(0);
    limits4->SetStats(0);
    
    //limits4->SetMarkerStyle(7);
    limits1->SetMarkerColor(1);
    //limits1->SetMarkerStyle(7);
    limits2->SetMarkerColor(2);
    //limits2->SetMarkerStyle(7);
    limits3->SetMarkerColor(3);
    //limits3->SetMarkerStyle(7);
    limits4->SetMarkerColor(9);
    //hs.Add(limits4);
    hs.Add(limits3);
    hs.Add(limits2);
    hs.Add(limits1);
	
	TCanvas * c2=new TCanvas("c2","",800,600);
	limits4->Draw("colz");	
    //hs.Draw("nostack");	
    c2->SaveAs((string("pdf_")+string(name)+string(".png")).c_str());
	
	TCanvas * c3=new TCanvas("c3","",800,600);
	limits1->Draw("colz");	
    //hs.Draw("nostack");	
    c3->SaveAs((string("pdf1_")+string(name)+string(".png")).c_str());
	
	TCanvas * c4=new TCanvas("c4","",800,600);
	limits2->Draw("colz");	
    //hs.Draw("nostack");	
    c4->SaveAs((string("pdf2_")+string(name)+string(".png")).c_str());


	TCanvas * c5=new TCanvas("c5","",800,600);
	limits3->Draw("colz");	
    //hs.Draw("nostack");	
    c5->SaveAs((string("pdf3_")+string(name)+string(".png")).c_str());

	TCanvas * c6=new TCanvas("c6","",800,600);
	limits5->Draw("colz");	
    //hs.Draw("nostack");	
    c6->SaveAs((string("pdf123_")+string(name)+string(".png")).c_str());
    
    TCanvas * c7=new TCanvas("c7","",800,600);
	hs.Draw("nostack");	
    c7->SaveAs((string("pdfall_")+string(name)+string(".png")).c_str());
	*/
	//cout<<"Lmax "<<lmax<<" "<<TMath::Prob(-2*log(lmax),m->size())<<" GaussMax "<<gaussmax<<endl;
	
	//cout<<"tanbmax "<<tanbmax<<" McHmax "<<mcHmax<<endl;
	cout<<"gL "<<gL<<" gQ "<<gQ<<endl;
	cout<<"lup "<<lup<<" qup "<<qup<<endl;
	
	cout<<"Btotaunu "<<m->BtotaunuR.subs(lst(m->tanb==tanbmax,m->McH==mcHmax))<<" exp "<<1.64/0.79<<" +/- "<<0.23*1.64/0.79<<endl;
	cout<<"BtoDtaunu "<<m->BtoDtaunuR.subs(lst(m->tanb==tanbmax,m->McH==mcHmax))<<" exp "<<440.0/296<<" +/- "<<1.4*58.0/296<<endl;
	cout<<"BtoD2taunu "<<m->BtoD2taunuR.subs(lst(m->tanb==tanbmax,m->McH==mcHmax))<<" exp "<<332.0/252<<" +/- "<<1.4*24.0/252<<endl;
	

	delete m;
}

	return 0;
	
}

