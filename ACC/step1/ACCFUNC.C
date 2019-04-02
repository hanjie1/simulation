// This is to simulate the x, Q2 shape of D2 and He3 kin1
// good electrons events
#define nTH 400
#define nNU 400
#define THmin 15.9
#define THmax 19.5
#define NUmin 7.36
#define NUmax 7.62
#define E0 10.589
#define Mp 0.938272 
#define nBin 25
void ACCFUNC()
{
     bool OUTINP=false;
     bool OUTYield=false;

     Double_t xbin[nBin];
     Double_t dBin=0.03;
     for(int ii=0;ii<nBin;ii++){
         xbin[ii]=0.15+ii*dBin;
     }

     Double_t XD2[nBin]={0.0},Q2_D2[nBin]={0.0},D2_Yield[nBin]={0.0};    
     Double_t XHe3[nBin]={0.0},Q2_He3[nBin]={0.0},He3_Yield[nBin]={0.0};    
     Double_t D2_err[nBin]={0.0},He3_err[nBin]={0.0};

     Double_t D2_pho=0.14215;  //g/cm^2
     Double_t He3_pho=0.0533752;  
     Double_t NA=TMath::Na();
 
     Double_t theta[nTH]={0.0},nu[nNU]={0.0};
     Double_t dTH=(THmax-THmin)/(nTH*1.0);
     for(int ii=0;ii<nTH;ii++)
	theta[ii]=THmin+ii*dTH;	

     Double_t dNU=(NUmax-NUmin)/(nNU*1.0);
     for(int ii=0;ii<nNU;ii++)
	nu[ii]=NUmin+ii*dNU;	

     Double_t wTH[nTH]={0.0},wNU[nNU]={0.0};
     for(int ii=0;ii<nTH;ii++){
	if(theta[ii]<16.5) wTH[ii]=5.0/3.0*theta[ii]-26.5;
	if(theta[ii]>=16.5 && theta[ii]<18.9) wTH[ii]=1.0;
	if(theta[ii]>=18.9) wTH[ii]=-5.0/3.0*theta[ii]+32.5;
     }
     TGraph *gTHACC=new TGraph(nTH,theta,wTH);

     for(int ii=0;ii<nNU;ii++){
	if(nu[ii]<7.4) wNU[ii]=25.0*nu[ii]-184;
	if(nu[ii]>=7.4 && nu[ii]<7.58) wNU[ii]=1.0;
	if(nu[ii]>=7.58) wNU[ii]=-25.0*nu[ii]+190.5;
     }
     TGraph *gNUACC=new TGraph(nNU,nu,wNU);

     TGraph2D *gXQ2=new TGraph2D(nNU*nTH);
     gXQ2->SetName("gXQ2");
     TGraph2D *gD2_XS=new TGraph2D(nNU*nTH);
     gD2_XS->SetName("gD2_XS");
     TGraph2D *gHe3_XS=new TGraph2D(nNU*nTH);
     gHe3_XS->SetName("gHe3_XS");

     TGraph *gR=new TGraph(nNU*nTH);
     gR->SetName("gR");

     int nn=0;

     TString outfile="ACCFUNC.inp";
     ofstream file;
     if(OUTINP){
       file.open(outfile);
       file<<"Marathon"<<endl;
       file<<outfile.Data()<<endl;
       file<<endl;
       file<<endl;
       file<<endl;
       file<<"E    Ep    theta"<<endl;
     }

     ifstream infile1;
     infile1.open("XStable/D2_xs.out");
     Ssiz_t from1=0;
     TString content1,tmp1;
     tmp1.ReadLine(infile1);

     ifstream infile2;
     infile2.open("XStable/He3_xs.out");
     Ssiz_t from2=0;
     TString content2,tmp2;
     tmp2.ReadLine(infile2);

     for(int ii=0;ii<nTH;ii++){
	for(int jj=0;jj<nNU;jj++){
	   Double_t tmp_weight=wNU[ii]*wTH[jj];
	   Double_t tmp_Q2=4.0*E0*(E0-nu[jj])*sin(theta[ii]*TMath::Pi()/(2.0*180))*sin(theta[ii]*TMath::Pi()/(2.0*180));
	   Double_t tmp_x=tmp_Q2/(2.0*Mp*nu[jj]);
	   gXQ2->SetPoint(nn,tmp_x,tmp_Q2,tmp_weight);
	   
	   tmp1.ReadLine(infile1);
	   tmp1.Tokenize(content1,from1," ");
	   Double_t in_x=atof(content1.Data());
	   if(abs(in_x-tmp_x)>0.0001){cout<<"Somethin wrong here for D2 x !!!  "<<in_x<<"  "<<tmp_x<<endl;break;}
	   tmp1.Tokenize(content1,from1," ");
	   Double_t in_Q2=atof(content1.Data());
	   if(abs(in_Q2-tmp_Q2)>0.0001){cout<<"Somethin wrong here for D2 Q2 !!!  "<<in_Q2<<"  "<<tmp_Q2<<endl;break;}
	   tmp1.Tokenize(content1,from1," ");
	   tmp1.Tokenize(content1,from1," ");
	   tmp1.Tokenize(content1,from1," ");
	   tmp1.Tokenize(content1,from1," ");
	   Double_t D2_rad=atof(content1.Data());
	   from1=0;
	
	   tmp2.ReadLine(infile2);
	   tmp2.Tokenize(content2,from2," ");
	   in_x=atof(content2.Data());
	   if(abs(in_x-tmp_x)>0.0001){cout<<"Somethin wrong here for He3 x !!!  "<<in_x<<"  "<<tmp_x<<endl;break;}
	   tmp2.Tokenize(content2,from2," ");
	   in_Q2=atof(content2.Data());
	   if(abs(in_Q2-tmp_Q2)>0.0001){cout<<"Somethin wrong here for He3 Q2 !!!  "<<in_Q2<<"  "<<tmp_Q2<<endl;break;}
	   tmp2.Tokenize(content2,from2," ");
	   tmp2.Tokenize(content2,from2," ");
	   tmp2.Tokenize(content2,from2," ");
	   tmp2.Tokenize(content2,from2," ");
	   Double_t He3_rad=atof(content2.Data());
	   from2=0;

	   for(int kk=0;kk<nBin;kk++){
	      if((tmp_x-xbin[kk])<dBin && (tmp_x-xbin[kk])>=0){
                 Double_t tmp_ND2=1.0/20.0;
                 Double_t tmp_NHe3=1.0/20.0;
		 XD2[kk]=XD2[kk]+tmp_x*tmp_weight*D2_rad*tmp_ND2; 
		 XHe3[kk]=XHe3[kk]+tmp_x*tmp_weight*He3_rad*tmp_NHe3; 
		 Q2_D2[kk]=Q2_D2[kk]+tmp_Q2*tmp_weight*D2_rad*tmp_ND2; 
		 Q2_He3[kk]=Q2_He3[kk]+tmp_Q2*tmp_weight*He3_rad*tmp_NHe3; 
		 D2_Yield[kk]+=tmp_weight*D2_rad*tmp_ND2;
		 He3_Yield[kk]+=tmp_weight*He3_rad*tmp_NHe3;
		 break;
	      }
	   }
	   
	   gD2_XS->SetPoint(nn,tmp_x,tmp_Q2,tmp_weight*D2_rad);
	   gHe3_XS->SetPoint(nn,tmp_x,tmp_Q2,tmp_weight*He3_rad);
	   gR->SetPoint(nn,tmp_x,He3_rad/D2_rad);
	   
	   nn++;
           if(OUTINP)file<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<E0-nu[jj]<<" "<<theta[ii]<<endl;
	}
     }
     infile1.close();
     infile2.close();
     if(OUTINP)file.close();

     TGraphErrors *gRatio=new TGraphErrors();    
     TGraph *gRatio_data=new TGraph();    
     TGraph *gRatio_diff=new TGraph();    

     Double_t He3_data[5]={32.1116028,26.5418167,19.4054222,14.7040148,12.2960033};
     Double_t D2_data[5]={20.2706814,16.7100353,12.1738443,9.19446468,7.67093134};
     int mm=0;
     ofstream outfile1;
     if(OUTYield){
       outfile1.open("Sim_Yield.dat");
       outfile1<<"x_D2    x_He3   Q2_D2    Q2_He3   D2_Yield   He3_Yield   ratio"<<endl;
     }
     for(int ii=0;ii<nBin;ii++){
	 if(XD2[ii]==0)continue;
	 Double_t xavg_D2=XD2[ii]/D2_Yield[ii];
	 Double_t xavg_He3=XHe3[ii]/He3_Yield[ii];
	 Double_t Q2avg_D2=Q2_D2[ii]/D2_Yield[ii];
	 Double_t Q2avg_He3=Q2_He3[ii]/He3_Yield[ii];
	 Double_t ratio=He3_Yield[ii]/D2_Yield[ii];
	 D2_err[ii]=ratio*sqrt(1.0/D2_Yield[ii]+1.0/He3_Yield[ii]);
 	 gRatio->SetPoint(ii,xavg_D2,ratio);
 	 gRatio->SetPointError(ii,0.0,D2_err[ii]);
         gRatio_data->SetPoint(mm,xavg_D2,He3_data[mm]/D2_data[mm]);
         gRatio_diff->SetPoint(mm,xavg_D2,(He3_data[mm]/D2_data[mm]-ratio)/ratio);
	 mm++;
	 if(OUTYield)outfile1<<xavg_D2<<"  "<<xavg_He3<<"  "<<Q2avg_D2<<"  "<<Q2avg_He3<<"  "<<D2_Yield[ii]<<"  "<<He3_Yield[ii]<<"  "<<ratio<<endl;
     }
     if(OUTYield)outfile1.close();
 
     TCanvas *c1=new TCanvas("c1","c1",1500,1500);
     gD2_XS->Draw("surf1");
     gD2_XS->SetTitle("MC Q2 vs x distribution");

     TCanvas *c2=new TCanvas("c2","c2",1500,1500);
     gHe3_XS->Draw("surf1");

     TCanvas *c3=new TCanvas("c3","c3",1500,1500);
     gXQ2->Draw("surf1");

     TCanvas *c4=new TCanvas("c4","c4",1500,1500);
     c4->Divide(2,1);
     c4->cd(1);
     TMultiGraph *mg1=new TMultiGraph();
     gRatio->SetMarkerStyle(8);
     gRatio->SetMarkerColor(4);
     gRatio_data->SetMarkerStyle(8);
     gRatio_data->SetMarkerColor(2);
     mg1->Add(gRatio);
     mg1->Add(gRatio_data);
     mg1->Draw("AP");
     mg1->SetTitle("Yield;x;");

     auto leg1=new TLegend(0.7,0.6,0.811,0.811);
     leg1->AddEntry(gRatio,"MC","P");
     leg1->AddEntry(gRatio_data,"rad_model","P");
     leg1->Draw();

     c4->cd(2);
     gRatio_diff->SetMarkerStyle(8);
     gRatio_diff->Draw("AP");
     gRatio_diff->SetTitle("(model-MC)/MC;x");

     TCanvas *c5=new TCanvas("c5","c5",1500,1500);
     c5->Divide(2,1);
     c5->cd(1);
     gTHACC->SetMarkerStyle(8);
     gTHACC->Draw("AP");
     gTHACC->SetTitle("theta ACC;theta");
     c5->cd(2);
     gNUACC->SetMarkerStyle(8);
     gNUACC->Draw("AP");
     gNUACC->SetTitle("nu ACC;nu");

}
