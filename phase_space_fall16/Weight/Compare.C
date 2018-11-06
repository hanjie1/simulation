#include "SetCut.h"
void Compare(){
     TFile *f1 = new TFile("Rootfiles/carbon_foil.root");
     TFile *f2 = new TFile("/home/hanjie/work/MARATHON/Rootfiles/tritium_1207.root");

     TTree *t1 =(TTree*)f1->Get("h1");
     TTree *t2 =(TTree*)f2->Get("T");

     TH1F *Dtg_th=new TH1F("Dtg_th","DATA theta at target",100,-0.06,0.06);
     TH1F *Mtg_th=new TH1F("Mtg_th","MC theta at target",100,-0.06,0.06);

     Float_t y_tg_rec,dp_rec,ph_tg_rec,th_tg_rec;
     Float_t born,crad,weight;
     t1->SetBranchAddress("hsytar",&y_tg_rec);
     t1->SetBranchAddress("hsdelta",&dp_rec);
     t1->SetBranchAddress("hsyptar",&ph_tg_rec);
     t1->SetBranchAddress("hsxptar",&th_tg_rec);
     t1->SetBranchAddress("xs_born",&born);
     t1->SetBranchAddress("xs_rad",&crad);
     t1->SetBranchAddress("weight",&weight);

     Int_t nentries=t1->GetEntries();
     for(int ii=0;ii<nentries;ii++){
	 t1->GetEntry(ii);
         if(abs(ph_tg_rec)<0.04 && abs(th_tg_rec)<0.06 && abs(dp_rec/100.0)<0.04){
	    Mtg_th->Fill(th_tg_rec,born*weight);
         }
     }

     TCanvas *c1=new TCanvas("c1","c1",1500,1500);
     Mtg_th->Draw();
     Mtg_th->SetLineColor(4);
     Double_t inte1=Mtg_th->Integral();
     cout<<inte1<<endl;
     t2->Draw("L.tr.tg_th>>Dtg_th",ACC+CK+Ep+trigger2+TRK+beta,"same"); 
     Double_t inte2=Dtg_th->Integral();
     //Dtg_th->Scale(inte1/inte2);
     Dtg_th->SetLineColor(2);


}
