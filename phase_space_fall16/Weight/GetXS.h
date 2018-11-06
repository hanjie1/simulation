#include <TMath.h>
const Double_t pi=TMath::Pi();
const int nTh=51;
const int nEp=51;
const Float_t delth=0.1;
const Float_t delep=0.005;
TString filename="XStable/Bodek/Carbon_kin1_xs.out";

class Get_XS
{
      public:
             Get_XS(){
                        //do nothing
             }
             virtual ~Get_XS(){
                        // delete aXS;
             }
       
      void LoadXStable(){
           ifstream infile;
           infile.open(filename);
           if(infile.is_open())cout<<"Open XS table"<<endl;
	   else {cout<<"Couldn't find XS table !!"<<endl;return;}

           int nn=0,mm=0;
     	   Ssiz_t from=0;
   	   TString content,tmp;
	   tmp.ReadLine(infile);
	   while(tmp.ReadLine(infile)){
                 tmp.Tokenize(content,from," ");
                 tmp.Tokenize(content,from," ");
         	 tmp.Tokenize(content,from," ");
           	 Theta_Table[mm]=atof(content.Data());
           	 tmp.Tokenize(content,from," ");
          	 Ep_Table[nn]=atof(content.Data());
        	 tmp.Tokenize(content,from," ");
          	 XSborn_Table[nn][mm]=atof(content.Data());
          	 tmp.Tokenize(content,from," ");
         	 XSrad_Table[nn][mm]=atof(content.Data());
          	 if((mm+1)%51==0){
		     mm=0;
		     nn++;
		 }
		 else{
		     mm++;
		 }
           	 from=0;
           }
     	   infile.close();
      }

      void SearchXS(Float_t Ep, Float_t aTheta, Float_t born, Float_t crad){
           Double_t th=aTheta*180.0/pi;
           int nfound=0;
	   int bbp=0,bbTh=0;
           born=0.0,crad=0.0;
	   Float_t aEp=Ep/1000.0;
//      theta
//	^         a3_______a4
//	|	   |___|___|     the XS table should be in ep increase and theta increase
//	|----->Ep  |___|___|
//                 a1     a2  
           for(int ii=0;ii<nEp;ii++){
               Double_t dEp=aEp-Ep_Table[ii];
               if(dEp<=delep && dEp>=0){
		  bbp=ii;
                  for(int jj=0;jj<nTh;jj++){
       	             Double_t dTheta=aTheta-Theta_Table[jj];
                     if(dTheta<=delth && dTheta>=0){
                        nfound=1;
		     }
		     if(nfound==1){bbTh=jj;break;}
	          }
                  if(nfound==0){cout<<"Theta is out of XS table range !!!!"; return;}
               }  
               if(nfound==1)break;
	   }     
           if(nfound==0){cout<<"Ep: "<<aEp<<" is out of XS table range !!!!"; return;}

           Float_t a1,a2,a3,a4;
           a1=XSborn_Table[bbp][bbTh];
	   a2=XSborn_Table[bbp+1][bbTh];
	   a3=XSborn_Table[bbp][bbTh+1];
	   a4=XSborn_Table[bbp+1][bbTh+1];

           Float_t tmp1=0.0,tmp2=0.0;
	   tmp1=a1+(a2-a1)/(Ep_Table[bbp+1]-Ep_Table[bbp])*(aEp-Ep_Table[bbp]);
	   tmp2=a3+(a4-a3)/(Ep_Table[bbp+1]-Ep_Table[bbp])*(aEp-Ep_Table[bbp]);

	   born=tmp1+(tmp2-tmp1)/(Theta_Table[bbTh+1]-Theta_Table[bbTh])*(aTheta-Theta_Table[bbTh]);

	   a1=XSrad_Table[bbp][bbTh];
	   a2=XSrad_Table[bbp+1][bbTh];
	   a3=XSrad_Table[bbp][bbTh+1];
	   a4=XSrad_Table[bbp+1][bbTh+1];

	   tmp1=a1+(a2-a1)/(Ep_Table[bbp+1]-Ep_Table[bbp])*(aEp-Ep_Table[bbp]);
	   tmp2=a3+(a4-a3)/(Ep_Table[bbp+1]-Ep_Table[bbp])*(aEp-Ep_Table[bbp]);

	   crad=tmp1+(tmp2-tmp1)/(Theta_Table[bbTh+1]-Theta_Table[bbTh])*(aTheta-Theta_Table[bbTh]);
           return;
      }          



      private:
	Double_t Ep_Table[nEp];
	Double_t Theta_Table[nTh];
	Double_t XSborn_Table[nEp][nTh];
	Double_t XSrad_Table[nEp][nTh];


};
