  Double_t th_max=0.06;
  Double_t th_min=-0.06;
  Double_t ph_max=0.03;
  Double_t ph_min=-0.03;
  Double_t dp_max=0.040;
  Double_t dp_min=-0.040;
  Double_t vz_max=0.10;
  Double_t vz_min=-0.09;

  TCut beta = "L.tr.beta>0";
  TCut Ep = "(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>0.75";
  TCut ACC = Form("L.tr.tg_th>%f && L.tr.tg_th<%f && L.tr.tg_ph>%f && L.tr.tg_ph<%f && L.tr.tg_dp>%f && L.tr.tg_dp<%f",th_min,th_max,ph_min,ph_max,dp_min,dp_max);
  TCut VZ = Form("L.tr.vz>%f && L.tr.vz<%f",vz_min,vz_max);
  TCut CK = "L.cer.asum_c>2000";
  TCut trigger2 = "(DL.evtypebits>>2)&1";
  TCut trigger1 = "(DL.evtypebits>>1)&1";
  TCut TRK = "L.tr.n==1";
  TCut MC_ACC=Form("hsxptar>%f && hsxptar<%f && hsyptar>%f && hsyptar<%f && hsdelta/100.0>%f && hsdelta/100.0<%f",th_min,th_max,ph_min,ph_max,dp_min,dp_max);
