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

  YAML::Node node = YAML::LoadFile("config_splitting.yaml");

  int hqpdg = node["hqpdg"].as<int>();
  int maxnevents = node["maxneventsperjob"].as<int>();
  int tune = node["tune"].as<int>();
  int beamidA = node["beamidA"].as<int>();
  int beamidB = node["beamidB"].as<int>();
  float eCM = node["eCM"].as<float>();
  float pTHatMin = node["pTHatMin"].as<float>();
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
  pythia.readString(Form("PhaseSpace:pTHatMin = %f", pTHatMin));
  if (FSR) pythia.readString("PartonLevel:FSR = on");
  else pythia.readString("PartonLevel:FSR = off");
  if (ISR) pythia.readString("PartonLevel:ISR = on");
  else pythia.readString("PartonLevel:ISR = off");

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
  //scatter plot kinematics
  TH2F *hPt1Pt2 = new TH2F("hPt1Pt2", ";p^{c2}_{T} (GeV/c);p^{c1}_{T} (GeV/c)",
		           200, 0., 200., 200., 0, 200.);
  TH2F *hPt1vsPtg = new TH2F("hPt1vsPtg", ";p^{g}_{T} (GeV/c);p^{c1}_{T} (GeV/c)",
		             200, 0., 200., 200., 0, 400.);
  TH1F *hPxConservat = new TH1F("hPxConservat", ";conservation gluon px;Entries",
                                1000., 0, 3.);
  //plots vs gluon E
  TH2F *hdeltaRvsGluonE = new TH2F("hdeltaRvsGluonE", ";E_{g} (GeV);#Delta R (c#bar{c})",
	                           1000, 0., 100., 100., 0, 3.);
  TH2F *hInvmassvsGluonE = new TH2F("hInvmassvsGluonE", ";E_{g} (GeV);mass(c#bar{c}) (GeV/c^{2})",
                                    200, 0., 200., 200., 0, 200.);
  TH2F *hGformtimevsGluonE = new TH2F("hGformtimevsGluonE", ";E_{g} (GeV); formation time(fm/c)",
		                      200, 0., 200., 200., 0, 20.);
  //study vs formation time
  TH2F *hGformtimevsDeltaR = new TH2F("hGformtimevsDeltaR", ";#Delta R; formation time(fm/c)",
		                      300, 0., 3., 200., 0, 20.);
  TH2F *hPtcharmvsGformtime = new TH2F("hPtcharmvsGformtime", ";formation time (fm/c); p^{c1}_{T} (GeV)",
		                       100, 0., 10., 500., 0, 100.);
  //single var distribution
  TH1F *hInvmass = new TH1F("hInvmass", ";mass (c#bar{c}) (GeV/c^{2};Entries",
		            100., 0, 20.);
  TH1F *hdeltaR = new TH1F("hdeltaR", ";#Delta R (c#bar{c});Entries", 1000., 0, 3.);
  TH1F *hetaG = new TH1F("hetaG", ";#eta(g);Entries", 1000., -5, 5.);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < maxnevents; ++iEvent) {
    pythia.next();
    for (int i = 0; i < pythia.event.size(); ++i) {
      // FIXME this is a consistency check on the calculation of the particle
      // energy I dont yet know what it does
      double eCalc = pythia.event[i].eCalc();
      if (abs(eCalc / pythia.event[i].e() - 1.) > 1e-3)
        cout << " e mismatch, i = " << i
             << " e_nominal = " << pythia.event[i].e()
             << " e-from-p = " << eCalc << " m-from-e "
             << pythia.event[i].mCalc() << "\n";
      if (pythia.event[i].pT() < 0 || pythia.event[i].pT() > 1.e+5)
        continue;
      // Selecting heavy quarks. Only quarks, not antiquark (4 for charm, 5 for
      // beauty)
      if (pythia.event[i].id() == hqpdg) {
        // checking the mothers of the heavy quarks. Remember there are two
        // mothers for each particle in Pythia
        int mothercharmFirstIndex = pythia.event[i].mother1();
        int mothercharmSecIndex = pythia.event[i].mother2();
        int pdgFirstMother = pythia.event[mothercharmFirstIndex].id();
        if (pdgFirstMother != 21)
          continue;
        int gluonDaugthFirstIndex =
            pythia.event.daughterList(mothercharmFirstIndex)[0];
        int gluonDaugthSecIndex =
            pythia.event.daughterList(mothercharmFirstIndex)[1];
        int pdggluonDaugthFirstIndex = pythia.event[gluonDaugthFirstIndex].id();
        int pdggluonDaugthSecIndex = pythia.event[gluonDaugthSecIndex].id();
        if (!((pdggluonDaugthFirstIndex == -hqpdg &&
               pdggluonDaugthSecIndex == hqpdg) ||
              (pdggluonDaugthFirstIndex == hqpdg &&
               pdggluonDaugthSecIndex == -hqpdg)))
          continue;
        int ndaughters = gluonDaugthSecIndex - gluonDaugthFirstIndex + 1;
        if (ndaughters != 2) {
          std::cout << "ATTENTION: the gluon has more than two daugthers"
                    << ndaughters << std::endl;
          continue;
        }
        // FIXME selection on gluon eta now hardcoded
        bool matchingparton =
            (i == gluonDaugthFirstIndex) || (i == gluonDaugthSecIndex);
        if (!matchingparton)
          continue;
        TLorentzVector vg, v1, v2;
        vg.SetPxPyPzE(pythia.event[mothercharmFirstIndex].px(),
		      pythia.event[mothercharmFirstIndex].py(),
		      pythia.event[mothercharmFirstIndex].pz(),
		      pythia.event[mothercharmFirstIndex].e());
        v1.SetPxPyPzE(pythia.event[gluonDaugthFirstIndex].px(),
		      pythia.event[gluonDaugthFirstIndex].py(),
		      pythia.event[gluonDaugthFirstIndex].pz(),
		      pythia.event[gluonDaugthFirstIndex].e());
        v2.SetPxPyPzE(pythia.event[gluonDaugthSecIndex].px(),
		      pythia.event[gluonDaugthSecIndex].py(),
		      pythia.event[gluonDaugthSecIndex].pz(),
		      pythia.event[gluonDaugthSecIndex].e());
        TLorentzVector vgfromcharms = v1 + v2;
	TVector3 vg3 = vg.Vect();
        TVector3 v1g = v1.Vect();
        TVector3 v2g = v2.Vect();
	TVector3 vgfromcharms3 = vgfromcharms.Vect();
        
        if (vgfromcharms.E() < minEgluon ||
            vgfromcharms.E() > maxEgluon ||
            std::abs(pythia.event[mothercharmFirstIndex].eta()) > maxabsEtaGluon)
          continue;
        // pythia.event.list();
        // std::cout<<"charm index"<<gluonDaugthFirstIndex<<std::endl;
        // std::cout<<"charm index"<<gluonDaugthSecIndex<<std::endl;
        // return 0;
        if (pythia.event[gluonDaugthFirstIndex].status() != -51 ||
            pythia.event[gluonDaugthFirstIndex].status() != -51) {
          std::cout << "ERROR: Status selected topologies"
                    << pythia.event[gluonDaugthFirstIndex].status()
                    << std::endl;
        }
	//FIXME: from here below calculations are made
	double pt1 = v1.Pt();
	double pt2 = v2.Pt();
	double ptg = vgfromcharms.Pt();
	double eta1 = v1.Eta();
	double eta2 = v2.Eta();
	double phi1 = v1.Phi();
	double phi2 = v2.Phi();
        double r =
            sqrt((eta2 - eta1) * (eta2 - eta1) + (phi2 - phi1) * (phi2 - phi1));
        double pxconserv = (v1.Px() + v2.Px() - vg.Px())/vg.Px();
        double Gformtime_thr = 0.2 * pythia.event[gluonDaugthFirstIndex].eCalc() /
                                (2 * 1.5*1.5);
        double Q2gluonfromtwocharms = vgfromcharms.M2();
        double Qgluonfromtwocharms = vgfromcharms.M();
        //std::cout<<"diff % E from gluon and from sum of charms=" 
	//	 <<(vgfromcharms.E()-vg.E())/vg.E()<<","<<vgfromcharms.E()<<","<<vg.E()<<std::endl;
        double Egluonfromtwocharms = vgfromcharms.E();
        double Gformtime =  0.2 *  Egluonfromtwocharms/Q2gluonfromtwocharms;
	
	//FIXME: from here below only histogram filling
        hPt1Pt2->Fill(pt1, pt2);
        hPt1vsPtg->Fill(ptg, pt1);
        hPxConservat->Fill(pxconserv);
        hetaG->Fill(vgfromcharms3.Eta());
        hdeltaR->Fill(r);
        hdeltaRvsGluonE->Fill(Egluonfromtwocharms, r);
        hInvmassvsGluonE->Fill(Egluonfromtwocharms, Qgluonfromtwocharms);
        hInvmass->Fill(Qgluonfromtwocharms);
        hGformtimevsGluonE->Fill(Egluonfromtwocharms, Gformtime);
	hGformtimevsDeltaR->Fill(r, Gformtime);
        hPtcharmvsGformtime->Fill(Gformtime, pt1);
      }
      // End of event loop. Statistics. Histogram. Done.
    }
  }
  pythia.stat();
  printf("nAccepted %ld, nTried %ld\n", pythia.info.nAccepted(),
         pythia.info.nTried());
  printf("pythia.info.sigmaGen() %f\n", pythia.info.sigmaGen());
  hetaG->Write();
  hdeltaR->Write();
  hdeltaRvsGluonE->Write();
  hPt1vsPtg->Write();
  TProfile *pPt1vsPtg =
      (TProfile *)hPt1vsPtg->ProfileX("pPt1vsPtg");
  pPt1vsPtg->GetYaxis()->SetTitle("p^{c1}_{T} (GeV/c)");
  pPt1vsPtg->GetXaxis()->SetTitle("p^{g}_{T} (GeV/c)");
  pPt1vsPtg->Write();
  hPt1Pt2->Write();
  TProfile *pPt1Pt2 = (TProfile *)hPt1Pt2->ProfileX("pPt1Pt2");
  pPt1Pt2->GetYaxis()->SetTitle("p^{c1}_{T} (GeV/c)");
  pPt1Pt2->GetXaxis()->SetTitle("p^{c2}_{T} (GeV/c)");
  pPt1Pt2->Write();
  hInvmassvsGluonE->Write();
  TProfile *pInvmassvsGluonE = (TProfile *)hInvmassvsGluonE->ProfileX("pInvmassvsGluonE");
  pInvmassvsGluonE->GetYaxis()->SetTitle("m_{c#bar{c}} (GeV)");
  pInvmassvsGluonE->GetXaxis()->SetTitle("E_{g} (GeV)");
  hInvmass->Write();
  TProfile *pGformtimevsGluonE =
      (TProfile *)hGformtimevsGluonE->ProfileX("pGformtimevsGluonE");
  hGformtimevsGluonE->Write();
  pGformtimevsGluonE->GetYaxis()->SetTitle("Formation time (fm/c)");
  pGformtimevsGluonE->GetXaxis()->SetTitle("E_{g} (GeV)");
  pGformtimevsGluonE->Write();
  hPtcharmvsGformtime->Write();
  TProfile *pPtcharmvsGformtime =
      (TProfile *)hPtcharmvsGformtime->ProfileX("pPtcharmvsGformtime");
  pPtcharmvsGformtime->Write();
  hPxConservat->Write();
  hGformtimevsDeltaR->Write();
  TProfile *pGformtimevsDeltaR =
      (TProfile *)hGformtimevsDeltaR->ProfileX("pGformtimevsDeltaR");
  pGformtimevsDeltaR->GetYaxis()->SetTitle("Formation time (fm/c)");
  pGformtimevsDeltaR->GetXaxis()->SetTitle("#Delta R(c\bar{c})");

  fout->Write();
  return 0;
}
