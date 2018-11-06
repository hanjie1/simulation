#include <TMath.h>
TString targetfile="../infiles/carbon_foil.inp";
TString runfile="RunInfo/Run1207.dat";
const Double_t Na=TMath::Na();
const Double_t Qe=TMath::Qe();

class Get_Runinfo
{
      public:
             Get_Runinfo(){
                        //do nothing
             }
             virtual ~Get_Runinfo(){
                        // delete aXS;
             }

	  Float_t GetPhaseSpace(){
		  ifstream infile;
		  infile.open(targetfile);
		  if(infile.is_open())cout<<"Open target file"<<endl;
		  else {cout<<"Couldn't open target file !!"<<endl; return 0;}

		  Float_t ngen=0.0,gdp_lo=0.0,gdp_hi=0.0,gth_lo=0.0,gth_hi=0.0,gphi_lo=0.0,gphi_hi=0.0;
		  Float_t Ep_central=0.0;

                  int nn=0,mm=0;
       	          Ssiz_t from=0;
                  TString content,tmp;
	          tmp.ReadLine(infile);
		  while(tmp(0)=='!')
		     tmp.ReadLine(infile);

	          tmp.Tokenize(content,from," "); from=0;
		  ngen=atof(content.Data());
	          tmp.ReadLine(infile);
	          tmp.Tokenize(content,from," "); from=0;
		  Ep_central=atof(content.Data());
	          tmp.ReadLine(infile);
	          tmp.ReadLine(infile);
	          tmp.Tokenize(content,from," "); from=0;
		  gdp_lo=atof(content.Data());
	          tmp.ReadLine(infile);
	          tmp.Tokenize(content,from," "); from=0;
		  gdp_hi=atof(content.Data());
	          tmp.ReadLine(infile);
	          tmp.Tokenize(content,from," "); from=0;
		  gth_lo=atof(content.Data());
	          tmp.ReadLine(infile);
	          tmp.Tokenize(content,from," "); from=0;
		  gth_hi=atof(content.Data());
	          tmp.ReadLine(infile);
	          tmp.Tokenize(content,from," "); from=0;
		  gphi_lo=atof(content.Data());
	          tmp.ReadLine(infile);
	          tmp.Tokenize(content,from," "); from=0;
		  gphi_hi=atof(content.Data());

                  infile.close();
		 
		  Float_t phasespace=0.0;
		  Ep_central=Ep_central/1000.0;
		  phasespace=(gth_hi-gth_lo)/1000.0*(gphi_hi-gphi_lo)/1000.0*(gdp_hi-gdp_lo)/100.0*Ep_central;
		  phasespace=phasespace/ngen;
		  return phasespace;
	     }

             void GetLum(int run_number,Float_t &LT, Float_t &LTerr, Float_t &Lum){
		     TString table="LHRStest";
     		     TSQLServer* Server = TSQLServer::Connect("mysql://halladb/triton-work","triton-user","3He3Hdata");

     		     TString  query=Form("select * from %s where run_number=%d;",table.Data(),run_number);
                     TSQLResult* result=Server->Query(query.Data());
     		     TSQLRow *row;
                     int nrows = result->GetRowCount();
		     Float_t Charge,target_thickness;
        	     LT=0;
      	             LTerr=0;
		     Charge=0.0;
		     target_thickness=0;
                     if (nrows==0) {
        		cout<<"run "<<run_number<<" is NOT in "<<table<<endl;
        		Server->Close();
        		return;
     		     }

		     string target;
     		     query=Form("select target,livetime,LT_err,charge,target_thickness from %s where run_number=%d;",table.Data(),run_number);
     		     result=Server->Query(query.Data());
     		     row=result->Next();
    		     target=row->GetField(0);
    		     LT=atof(row->GetField(1));
     		     LTerr=atof(row->GetField(2));
    		     Charge=atof(row->GetField(3));
     		     target_thickness=atof(row->GetField(4));
		     Float_t massA=1.0;
		     if(target=="C")massA=12.0;
		     if(target=="He3")massA=3.0;
		     if(target=="H3")massA=3.0;
		     if(target=="D2")massA=2.0;
                     Server->Close();
                     Lum=target_thickness*Na/massA*Charge/(1.e33*Qe*1e6);  // 1/nb
		     return;
	     }	



};
