#include "GetXS.h"
#include "GetRuninfo.h"
void Weight()
{
     TFile *f1 = new TFile("Rootfiles/carbon_foil_nooffset_f1f217.root","update");
     assert(f1);

     TTree *t1 =(TTree*)f1->Get("h1");  
     assert(t1);

     Float_t x_fp,y_fp,th_fp,ph_fp;
     Float_t y_tg_gen,dp_gen,ph_tg_gen,th_tg_gen;
     Float_t y_tg_rec,dp_rec,ph_tg_rec,th_tg_rec;
     Float_t beam_x,beam_y;
     Float_t ok_spec,stopwhen;
     Float_t reactz_gen,reactz_rec,theta_gen,theta_rec,p_gen,p_rec,E_gen,E_rec,elossi,elossf;
     t1->SetBranchAddress("hsxfp",&x_fp);
     t1->SetBranchAddress("hsyfp",&y_fp);
     t1->SetBranchAddress("hsxpfp",&th_fp);
     t1->SetBranchAddress("hsypfp",&ph_fp);

     t1->SetBranchAddress("hsytari",&y_tg_gen);
     t1->SetBranchAddress("hsdeltai",&dp_gen);
     t1->SetBranchAddress("hsyptari",&ph_tg_gen);
     t1->SetBranchAddress("hsxptari",&th_tg_gen);

     t1->SetBranchAddress("hsytar",&y_tg_rec);
     t1->SetBranchAddress("hsdelta",&dp_rec);
     t1->SetBranchAddress("hsyptar",&ph_tg_rec);
     t1->SetBranchAddress("hsxptar",&th_tg_rec);

     t1->SetBranchAddress("frx",&beam_x);
     t1->SetBranchAddress("fry",&beam_y);
     t1->SetBranchAddress("ok_spec",&ok_spec);
     t1->SetBranchAddress("stopwhen",&stopwhen);

     t1->SetBranchAddress("z_init",&reactz_gen);
     t1->SetBranchAddress("z_recon",&reactz_rec);
     t1->SetBranchAddress("th_init",&theta_gen);
     t1->SetBranchAddress("th_recon",&theta_rec);
     t1->SetBranchAddress("p_init",&p_gen);
     t1->SetBranchAddress("p_recon",&p_rec);
     t1->SetBranchAddress("e_init",&E_gen);
     t1->SetBranchAddress("e_recon",&E_rec);
     t1->SetBranchAddress("elossi",&elossi);
     t1->SetBranchAddress("elossf",&elossf);

     Float_t xs_born,xs_rad,xbj,Q2,W2,nu,weight;
     TBranch *t1_born=t1->Branch("xs_born",&xs_born,"xs_born/F"); //born cross section
     TBranch *t1_rad=t1->Branch("xs_rad",&xs_rad,"xs_rad/F"); //radiative cross section
     TBranch *t1_xbj=t1->Branch("xbj",&xbj,"xbj/F"); //x bjokern
     TBranch *t1_nu=t1->Branch("nu",&nu,"nu/F"); //x bjokern
     TBranch *t1_q2=t1->Branch("Q2",&Q2,"Q2/F"); //Q2
     TBranch *t1_w2=t1->Branch("W2",&W2,"W2/F"); //W2
     TBranch *t1_weight=t1->Branch("weight",&weight,"weight/F"); //weight

     Get_XS *xs=new Get_XS();
     xs->LoadXStable();

     Float_t phasespace=0.0; 
     Get_Runinfo *Run=new Get_Runinfo();
     phasespace=Run->GetPhaseSpace();
     cout<<phasespace<<endl;

     Float_t LT,LTerr,Lum;
     int run_number=0;
     cout<<"Input run number: ";
     cin>>run_number;
     Run->GetLum(run_number,LT,LTerr,Lum);
     cout<<Lum<<"  "<<LT<<"  "<<LTerr<<"  "<<endl;

     Float_t E0=10.5901;
     Float_t Mp=0.938272;
     Int_t nentries=t1->GetEntries();
     for(int ii=0;ii<nentries;ii++){
         t1->GetEntry(ii);
         xs_born=0.0; xs_rad=0.0;	 
	 xbj=0.0;Q2=0.0;W2=0.0;nu=0.0;weight=0.0;
         if(ok_spec){
            xs->SearchXS(p_gen,theta_gen,xs_born,xs_rad); 	    
            if(xs_born==0||xs_rad==0){cout<<"Event "<<ii<<" doesn't have born XS"<<endl;continue;}
	    Q2=4.0*E0*E_rec/1000.0*sin(theta_rec/2.0)*sin(theta_rec/2.0);
	    nu=E0-E_rec/1000.0;
            xbj=Q2/(2.0*Mp*nu);
            W2=Mp*Mp+2*Mp*nu-Q2;        
	    weight=Lum*phasespace*LT;
	 }
//         cout<<xbj<<" "<<Q2<<" "<<nu<<" "<<W2<<" "<<weight<<" "<<xs_born<<" "<<xs_rad<<endl;
         t1_born->Fill();
         t1_rad->Fill();
         t1_xbj->Fill();
         t1_nu->Fill();
         t1_q2->Fill();
         t1_w2->Fill();
         t1_weight->Fill();

     }
     t1->Write();
     f1->Close();

}
