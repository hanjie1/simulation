#include "SetCut.h"
void Compare(){
     TFile *f1 = new TFile("Rootfiles/carbon_foil_nooffset_f1f217.root");
//     TFile *f3 = new TFile("Rootfiles/carbon_foil_nooffset_f1f217.root");
     TFile *f2 = new TFile("/home/hanjie/work/MARATHON/Rootfiles/tritium_1207.root");

     TTree *t1 =(TTree*)f1->Get("h1");
//     TTree *t3 =(TTree*)f3->Get("h1");
     TTree *t2 =(TTree*)f2->Get("T");

     TH1F *Mtg_th=new TH1F("Mtg_th","MC theta at target",100,-0.07,0.07);
     TH1F *Mtg_ph=new TH1F("Mtg_ph","MC phi at target",100,-0.05,0.05);
     TH1F *Mtg_dp=new TH1F("Mtg_dp","MC dp at target",100,-5,5);
     TH1F *Mtg_y=new TH1F("Mtg_y","MC target y",100,-2,2);
     TH1F *Mtg_z=new TH1F("Mtg_z","MC target z",100,-1.5,3);

     Float_t y_tg_rec,dp_rec,ph_tg_rec,th_tg_rec,z_recon;
     Float_t born,crad,weight;
     t1->SetBranchAddress("hsytar",&y_tg_rec);
     t1->SetBranchAddress("hsdelta",&dp_rec);
     t1->SetBranchAddress("hsyptar",&ph_tg_rec);
     t1->SetBranchAddress("hsxptar",&th_tg_rec);
     t1->SetBranchAddress("xs_born",&born);
     t1->SetBranchAddress("xs_rad",&crad);
     t1->SetBranchAddress("weight",&weight);
     t1->SetBranchAddress("z_recon",&z_recon);

     Int_t nentries=t1->GetEntries();
     for(int ii=0;ii<nentries;ii++){
	 t1->GetEntry(ii);
         if(abs(ph_tg_rec)<0.04 && abs(th_tg_rec)<0.06 && abs(dp_rec/100.0)<0.04){
	    Mtg_th->Fill(th_tg_rec,crad*weight);
	    Mtg_ph->Fill(ph_tg_rec,crad*weight);
	    Mtg_dp->Fill(dp_rec,crad*weight);
	    Mtg_y->Fill(y_tg_rec,crad*weight);
	    Mtg_z->Fill(z_recon,crad*weight);
         }
     }
/*
     TH1F *Mtg1_th=new TH1F("Mtg1_th","MC theta at target",100,-0.07,0.07);
     TH1F *Mtg1_ph=new TH1F("Mtg1_ph","MC phi at target",100,-0.05,0.05);
     TH1F *Mtg1_dp=new TH1F("Mtg1_dp","MC dp at target",100,-5,5);
     TH1F *Mtg1_y=new TH1F("Mtg1_y","MC target y",100,-3,3);

     Float_t y_tg1_rec,dp_rec1,ph_tg1_rec,th_tg1_rec;
     Float_t born1,crad1,weight1;
     t3->SetBranchAddress("hsytar",&y_tg1_rec);
     t3->SetBranchAddress("hsdelta",&dp_rec1);
     t3->SetBranchAddress("hsyptar",&ph_tg1_rec);
     t3->SetBranchAddress("hsxptar",&th_tg1_rec);
     t3->SetBranchAddress("xs_born",&born1);
     t3->SetBranchAddress("xs_rad",&crad1);
     t3->SetBranchAddress("weight",&weight1);

     nentries=t3->GetEntries();
     for(int ii=0;ii<nentries;ii++){
	 t3->GetEntry(ii);
         if(abs(ph_tg1_rec)<0.04 && abs(th_tg1_rec)<0.06 && abs(dp_rec1/100.0)<0.04){
	    Mtg1_th->Fill(th_tg1_rec,crad1*weight1);
	    Mtg1_ph->Fill(ph_tg1_rec,crad1*weight1);
	    Mtg1_dp->Fill(dp_rec1,crad1*weight1);
	    Mtg1_y->Fill(y_tg1_rec,crad1*weight1);
         }
     }
*/
     TH1F *Dtg_th=new TH1F("Dtg_th","Data theta at target",100,-0.07,0.07);
     TH1F *Dtg_ph=new TH1F("Dtg_ph","Data phi at target",100,-0.05,0.05);
     TH1F *Dtg_dp=new TH1F("Dtg_dp","Data dp at target",100,-5,5);
     TH1F *Dtg_y=new TH1F("Dtg_y","Data target y",100,-2,2);
     TH1F *Dtg_z=new TH1F("Dtg_z","Data target z",100,-1.5,3);

     TCanvas *c1=new TCanvas("c1","c1",1500,1500);
     Mtg_th->Draw();
     Mtg_th->SetLineColor(4);
//     Mtg1_th->Draw("SAME");
//     Mtg1_th->SetLineColor(1);
     t2->Draw("L.tr.tg_th>>Dtg_th",ACC+CK+Ep+trigger2+TRK+beta,"same"); 
     Dtg_th->SetLineColor(2);

     auto leg1 = new TLegend(0.7,0.7,0.85,0.85);
     leg1->AddEntry(Dtg_th,"Data","L");
     leg1->AddEntry(Mtg_th,"MC","L");
//     leg1->AddEntry(Mtg1_th,"f1f217","L");
     leg1->Draw();


     TCanvas *c2=new TCanvas("c2","c2",1500,1500);
     Mtg_ph->Draw();
     Mtg_ph->SetLineColor(4);
//     Mtg1_ph->Draw("same");
//     Mtg1_ph->SetLineColor(1);
     t2->Draw("L.tr.tg_ph>>Dtg_ph",ACC+CK+Ep+trigger2+TRK+beta,"same"); 
     Dtg_ph->SetLineColor(2);

     auto leg2 = new TLegend(0.7,0.7,0.85,0.85);
     leg2->AddEntry(Dtg_ph,"Data","L");
     leg2->AddEntry(Mtg_ph,"MC","L");
//     leg2->AddEntry(Mtg1_ph,"f1f217","L");
     leg2->Draw();

     TCanvas *c3=new TCanvas("c3","c3",1500,1500);
     Mtg_dp->Draw();
     Mtg_dp->SetLineColor(4);
//     Mtg1_dp->Draw("same");
//     Mtg1_dp->SetLineColor(1);
     t2->Draw("L.tr.tg_dp*100>>Dtg_dp",ACC+CK+Ep+trigger2+TRK+beta,"same"); 
     Dtg_dp->SetLineColor(2);

     auto leg3 = new TLegend(0.7,0.7,0.85,0.85);
     leg3->AddEntry(Dtg_dp,"Data","L");
     leg3->AddEntry(Mtg_dp,"MC","L");
//     leg3->AddEntry(Mtg1_dp,"f1f217","L");
     leg3->Draw();

     TCanvas *c4=new TCanvas("c4","c4",1500,1500);
     Mtg_y->Draw();
     Mtg_y->SetLineColor(4);
//     Mtg1_y->Draw("SAME");
//     Mtg1_y->SetLineColor(1);
     t2->Draw("L.tr.tg_y*100>>Dtg_y",ACC+CK+Ep+trigger2+TRK+beta,"same"); 
     Dtg_y->SetLineColor(2);

     auto leg4 = new TLegend(0.7,0.7,0.85,0.85);
     leg4->AddEntry(Dtg_y,"Data","L");
     leg4->AddEntry(Mtg_y,"MC","L");
//     leg4->AddEntry(Mtg1_y,"f1f217","L");
     leg4->Draw();

     TCanvas *c5=new TCanvas("c5","c5",1500,1500);
     Mtg_z->Draw();
     Mtg_z->SetLineColor(4);
//     Mtg1_z->Draw("SAME");
//     Mtg1_z->SetLineColor(1);
     t2->Draw("L.tr.vz*100>>Dtg_z",ACC+CK+Ep+trigger2+TRK+beta,"same"); 
     Dtg_z->SetLineColor(2);

     auto leg5 = new TLegend(0.7,0.7,0.85,0.85);
     leg5->AddEntry(Dtg_z,"Data","L");
     leg5->AddEntry(Mtg_z,"MC","L");
//     leg5->AddEntrz(Mtg1_z,"f1f217","L");
     leg5->Draw();
}
