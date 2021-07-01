// This is a simple macro to generate pythia predictions for the study
// of gluon splitting processes. The basic idea is to find g->ccbar
// topologies, identify the kinematic of the gluon and of the c/cbar
// and estimate the formation time of the gluon. In the future, a
// simple parametrization of the broadening process will be added.

#include "Pythia8/Pythia.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
#include "fastjet/PseudoJet.hh"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/ClusterSequenceArea.hh"
#include <cstdio> // needed for io
#include <ctime>
#include <iostream> // needed for io
#include <time.h>   /* time */
#include <valarray>
//#include <yaml.h>
// include <stdio.h>
// include <glib.h>
#include <yaml-cpp/yaml.h>

using namespace Pythia8;

int main(int argc, char *argv[]) {
  (void)argc;
  std::string mycase = argv[1];
  int cislo = -1; // unique number for each file
  cislo = atoi(argv[2]);
  // number of parallel jobs to be run. Be aware that each single file
  // will be normalized by this number to make sure that the merged output
  // file has the proper normalization!

  YAML::Node node = YAML::LoadFile("config_longrangecharm.yaml");

  int hqpdg = node["hqpdg"].as<int>();
  int maxnevents = node["maxneventsperjob"].as<int>();
  int tune = node["tune"].as<int>();
  int beamidA = node["beamidA"].as<int>();
  int beamidB = node["beamidB"].as<int>();
  float eCM = node["eCM"].as<float>();
  float minEgluon = node["minEgluon"].as<float>();
  float maxEgluon = node["maxEgluon"].as<float>();
  float maxabsEtaGluon = node["maxabsEtaGluon"].as<float>();
  const std::string pythiamode = node["pythiamode"].as<std::string>();
  const std::string outputfile = node["outputfile"].as<std::string>();
  const std::string extramode = node["extramode"].as<std::string>();
  bool ISR = node["ISR"].as<bool>();
  bool FSR = node["FSR"].as<bool>();

  // END OF CONFIGURATION

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString(Form("%s", pythiamode.data()));
  pythia.readString(Form("Main:numberOfEvents = %d", maxnevents));
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString(Form("Tune:pp = %d", tune));
  pythia.readString(Form("Beams:idA = %d", beamidA));
  pythia.readString(Form("Beams:idB = %d", beamidB));
  pythia.readString(Form("Beams:eCM = %f", eCM));
  if (FSR)
    pythia.readString("PartonLevel:FSR = on");
  else
    pythia.readString("PartonLevel:FSR = off");
  if (ISR)
    pythia.readString("PartonLevel:ISR = on");
  else
    pythia.readString("PartonLevel:ISR = off");

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d", cislo));
  if (extramode == "mode2") {
    std::cout << "Running with mode2" << std::endl;
    pythia.readString("ColourReconnection:mode = 1");
    pythia.readString("ColourReconnection:allowDoubleJunRem = off");
    pythia.readString("ColourReconnection:m0 = 0.3");
    pythia.readString("ColourReconnection:allowJunctions = on");
    pythia.readString("ColourReconnection:junctionCorrection = 1.20");
    pythia.readString("ColourReconnection:timeDilationMode = 2");
    pythia.readString("ColourReconnection:timeDilationPar = 0.18");
    pythia.readString("StringPT:sigma = 0.335");
    pythia.readString("StringZ:aLund = 0.36");
    pythia.readString("StringZ:bLund = 0.56");
    pythia.readString("StringFlav:probQQtoQ = 0.078");
    pythia.readString("StringFlav:ProbStoUD = 0.2");
    pythia.readString(
        "StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
    pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
    pythia.readString("BeamRemnants:remnantMode = 1");
    pythia.readString("BeamRemnants:saturation =5");
  }
  pythia.init();

  TFile *fout = new TFile(outputfile.data(), "recreate");

  fout->cd();
  // scatter plot kinematics
  TH1F *hDeltaPhi = new TH1F("hDeltaPhi", ";#Phi #eta;Entries", 80, -4., 4.);
  TH1F *hDeltaEta = new TH1F("hDeltaEta", ";#Delta #eta;Entries", 80, -4., 4.);
  TH1F *hDeltaPt = new TH1F("hDeltaPt", ";#Delta p_{T};Entries", 100, 0., 10.);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < maxnevents; ++iEvent) {
    pythia.next();
    vector<int> dplus, dminus;
    for (int i = 0; i < pythia.event.size(); ++i) {
      // FIXME this is a consistency check on the calculation of the particle
      // energy I dont yet know what it does
      
/*
      double eCalc = pythia.event[i].eCalc();
      if (abs(eCalc / pythia.event[i].e() - 1.) > 1e-3)
        cout << " e mismatch, i = " << i
             << " e_nominal = " << pythia.event[i].e()
             << " e-from-p = " << eCalc << " m-from-e "
             << pythia.event[i].mCalc() << "\n";
*/
       	if (pythia.event[i].pT() < 0 || pythia.event[i].pT() > 1.e+5)
        continue;
      // Selecting heavy quarks. Only quarks, not antiquark (4 for charm, 5 for
      // beauty)

      if (pythia.event[i].id() == 4) {std::cout<<"charm="<<pythia.event[i].pT()<<std::endl; dplus.push_back(i);}
      if (pythia.event[i].id() == -4) {std::cout<<"anticharm="<<pythia.event[i].pT()<<std::endl; dminus.push_back(i);}
    }
    int npar = dplus.size();
    int nantipar = dminus.size();
    if (npar == 1 && nantipar == 1) {
      hDeltaEta->Fill(pythia.event[dplus[0]].eta() - pythia.event[dminus[0]].eta());
      hDeltaPt->Fill(pythia.event[dplus[0]].pT() - pythia.event[dminus[0]].pT());
      hDeltaPhi->Fill(pythia.event[dplus[0]].phi() - pythia.event[dminus[0]].pT());
    }
  dplus.clear();
  dminus.clear();
  } // end of loop over events

  // Check whether pair(s) present.

  pythia.stat();
  printf("nAccepted %ld, nTried %ld\n", pythia.info.nAccepted(),
         pythia.info.nTried());
  printf("pythia.info.sigmaGen() %f\n", pythia.info.sigmaGen());
  hDeltaEta->Write();
  fout->Write();
  return 0;
}
