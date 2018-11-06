#include "GetXS.h"
#include "GetRuninfo.h"
void Weight()
{
     TFile *f1 = new TFile("../worksim/carbon_foil.root","update");
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

     Float_t xs_born,xs_rad,xbj,Q2,W2,weight;
     t1->Branch("xs_born",&xs_born,"xs_born/D"); //born cross section
     t1->Branch("xs_rad",&xs_rad,"xs_rad/D"); //radiative cross section
     t1->Branch("xbj",&xbj,"xbj/D"); //x bjokern
     t1->Branch("Q2",&Q2,"Q2/D"); //Q2
     t1->Branch("W2",&W2,"W2/D"); //W2
     t1->Branch("weight",&weight,"weight/D"); //weight

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

     Int_t nentries=t1->GetEntries();
     for(int ii=0;ii<100;ii++){
         xs_born=0.0; xs_rad=0.0;	 
         t1->GetEntry(ii);
         if(ok_spec){
	    cout<<theta_gen<<endl;
            xs->SearchXS(p_gen,theta_gen,xs_born,xs_rad); 	    
            // if(xs_born==0||xs_rad==0){cout<<"Event "<<ii<<" doesn't have born XS"<<endl;continue;}
            
	 }

     }



}
