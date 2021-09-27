void getRaa()
{

  TH1F *histCentral = new TH1F("RaaCentral", "Raa 0-10%; pT (GeV/c); Raa", 60, 0, 12);
  TH1F *histSemi = new TH1F("RaaSemi", "Raa 30-50%; pT (GeV/c); Raa", 60, 0, 12);

  //  read the tables
 
  double pTCentral[60] = {0};
  double RaaLowCentral[60] = {0};
  double RaaHighCentral[60] = {0};
  double pTSemi[60] = {0};
  double RaaLowSemi[60] = {0};
  double RaaHighSemi[60] = {0};
  
  double pT, RaaLow, RaaHigh = 0;
 
  ifstream readCentral;
  readCentral.open("../InputsTheory/LuukRaa/PbPb5.02TeV_Lc_RAA_010.txt");
  int i=0;
  while(readCentral >> pT)
  {
    readCentral >> RaaLow >> RaaHigh;
    //pTCentral[i] = pT;  
    //RaaLowCentral[i] = RaaLow;
    //RaaHighCentral[i] = RaaHigh;

    histCentral->SetBinContent(i+1, (RaaHigh+RaaLow)/2.0);
    histCentral->SetBinError(i+1, (RaaHigh-RaaLow)/2.0);

    i++;
  }

  ifstream readSemi;
  readSemi.open("../InputsTheory/LuukRaa/PbPb5.02TeV_Lc_RAA_3050.txt");
  int j=0;
  while(readSemi >> pT)
  {
    readSemi >> RaaLow >> RaaHigh;
    //pTSemi[j] = pT;  
    //RaaLowSemi[j] = RaaLow;
    //RaaHighSemi[j] = RaaHigh;

    histSemi->SetBinContent(j+1, (RaaHigh+RaaLow)/2.0);
    histSemi->SetBinError(j+1, (RaaHigh-RaaLow)/2.0);

    j++;
  }

  histCentral->Draw(); 
  histSemi->Draw("same"); 

  //  save into root files
  TFile *outCentral = new TFile("../InputsTheory/Lambda_c_ptdep_TAMU_PbPb5p02_cent010_absy0p5.root", "recreate");
  histCentral->Write();
  outCentral->Close();

  TFile *outSemi = new TFile("../InputsTheory/Lambda_c_ptdep_TAMU_PbPb5p02_cent3050_absy0p5.root", "recreate");
  histSemi->Write();
  outSemi->Close();

}
