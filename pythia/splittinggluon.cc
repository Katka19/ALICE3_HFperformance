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
  const std::string pythiamode = node["pythiamode"].as<std::string>();
  const std::string outputfile = node["outputfile"].as<std::string>();
  const std::string extramode = node["extramode"].as<std::string>();

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
  // pythia.readString("PartonLevel:FSR = off");
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

  // difference between gluon px and sum of charm px, normalized by gluon px
  TH1F *hPxConservat =
      new TH1F("hPxConservat", ";conservation gluon px in splitting;Entries",
               1000., 0, 3.);
  // deltaR between the c-cbar pair defined in terms of eta and phi
  TH1F *hdeltaR =
      new TH1F("hdeltaR", ";#Delta R (c#bar{c});Entries", 1000., 0, 3.);
  // distribution of gluon eta
  TH1F *hetaG = new TH1F("hetaG", ";#eta(g);Entries", 1000., -5, 5.);
  // deltaR between the c-cbar pair defined in terms of eta and phi as a
  // function of the gluon energy
  TH2F *hdeltaRvsGluonE =
      new TH2F("hdeltaRvsGluonE", ";E(g);#Delta R (c#bar{c})", 1000, 0., 100.,
               100., 0, 3.);
  // scatter plot pt charm vs pt anticharm
  TH2F *hPt1Pt2 = new TH2F("hPt1Pt2", ";p_{T}(c);p_{T}(cbar)", 200, 0., 200.,
                           200., 0, 200.);
  // scatter plot pt charm vs Egluon
  TH2F *hPtcharmvsGluonE = new TH2F("hPtcharmvsGluonE", ";E(g);p_{T}(c)", 200,
                                    0., 200., 200., 0, 200.);
  // scatter plot inv mass ccbar vs gluon energy
  TH2F *hInvmassvsGluonE = new TH2F("hInvmassvsGluonE", ";E(g);mass(c#bar{c})",
                                    200, 0., 200., 200., 0, 200.);
  // Inv. mass ccbar
  TH1F *hInvmass =
      new TH1F("hInvmass", ";mass (c#bar{c});Entries", 100., 0, 20.);
  // Gluon formation time vs Gluon energy
  TH2F *hGformtimevsGluonE =
      new TH2F("hGformtimevsGluonE", ";gluon energy; formation time(fm/c)", 200,
               0., 200., 200., 0, 200.);
  // perperdicular component of the 3 momentum c vs cbar vs gluon formation time
  TH2F *hPperpGformtime =
      new TH2F("hPperpGformtime", ";formation time (fm/c); q_{perp} (GeV)", 100,
               0., 10., 500., 0, 100.);
  // pT of the charm vs gluon formation time
  TH2F *hPtcharmvsGformtime = new TH2F(
      "hPtcharmvsGformtime", ";formation time (fm/c); p_{T} charm (GeV)", 100,
      0., 10., 500., 0, 100.);

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
        if ((mothercharmSecIndex - mothercharmFirstIndex) != 0 &&
            mothercharmSecIndex != 0) {
          std::cout << "ATTENTION: the charm quark has "
                    << mothercharmSecIndex - mothercharmFirstIndex + 1
                    << "mothers, it is not a splitting!" << std::endl;
          // pythia.event.list();
          // std::cout<<"charm index"<<gluonDaugthFirstIndex<<std::endl;
          // std::cout<<"charm index"<<gluonDaugthSecIndex<<std::endl;
          // return 0;
          // for (int ind=mothercharmSecIndex; ind<=mothercharmFirstIndex;
          // ind++){
          //  std::cout<<pythia.event[ind].id()<<std::endl;
          //}
          continue;
        }
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
        if (pythia.event[mothercharmFirstIndex].e() < minEgluon ||
            std::abs(pythia.event[mothercharmFirstIndex].eta()) > 1.)
          continue;
        // pythia.event.list();
        // std::cout<<"charm index"<<gluonDaugthFirstIndex<<std::endl;
        // std::cout<<"charm index"<<gluonDaugthSecIndex<<std::endl;
        // return 0;
        double pt1 = pythia.event[gluonDaugthFirstIndex].pT();
        double pt2 = pythia.event[gluonDaugthSecIndex].pT();
        double eta1 = pythia.event[gluonDaugthFirstIndex].eta();
        double eta2 = pythia.event[gluonDaugthSecIndex].eta();
        double phi1 = pythia.event[gluonDaugthFirstIndex].phi();
        double phi2 = pythia.event[gluonDaugthSecIndex].phi();
        double r =
            sqrt((eta2 - eta1) * (eta2 - eta1) + (phi2 - phi1) * (phi2 - phi1));
        double massccbar = (pythia.event[gluonDaugthFirstIndex].p() +
                            pythia.event[gluonDaugthSecIndex].p())
                               .mCalc();
        double pxconserv = pythia.event[gluonDaugthFirstIndex].px() +
                           pythia.event[gluonDaugthSecIndex].px() -
                           pythia.event[mothercharmFirstIndex].px();

        hPxConservat->Fill(pxconserv /
                           pythia.event[mothercharmFirstIndex].px());
        hetaG->Fill(pythia.event[mothercharmFirstIndex].eta());
        hdeltaR->Fill(r);
        hdeltaRvsGluonE->Fill(pythia.event[mothercharmFirstIndex].e(), r);
        hPtcharmvsGluonE->Fill(pythia.event[mothercharmFirstIndex].e(), pt1);
        hPt1Pt2->Fill(pt1, pt2);
        hInvmassvsGluonE->Fill(pythia.event[mothercharmFirstIndex].e(),
                               massccbar);
        hInvmass->Fill(massccbar);
        // hcut*c=197.3MeV fm
        // In NU, 1 = 0.2 GeV * fm
        // 5 GeV-1 = fm
        // GeV-1 = 0.2 fm
        // X GeV-1 = 0.2 X fm
        double Gformtime = 0.2 * pythia.event[gluonDaugthFirstIndex].eCalc() /
                           (2 * pythia.event[gluonDaugthFirstIndex].m2Calc());
        hGformtimevsGluonE->Fill(pythia.event[mothercharmFirstIndex].eCalc(),
                                 Gformtime);
        TLorentzVector vg(pythia.event[mothercharmFirstIndex].e(),
                          pythia.event[mothercharmFirstIndex].px(),
                          pythia.event[mothercharmFirstIndex].py(),
                          pythia.event[mothercharmFirstIndex].pz());
        TLorentzVector v1(pythia.event[gluonDaugthFirstIndex].e(),
                          pythia.event[gluonDaugthFirstIndex].px(),
                          pythia.event[gluonDaugthFirstIndex].py(),
                          pythia.event[gluonDaugthFirstIndex].pz());
        TLorentzVector v2(pythia.event[gluonDaugthSecIndex].e(),
                          pythia.event[gluonDaugthSecIndex].px(),
                          pythia.event[gluonDaugthSecIndex].py(),
                          pythia.event[gluonDaugthSecIndex].pz());
        TVector3 pg = vg.Vect();
        TVector3 p1 = v1.Vect();
        TVector3 p2 = v2.Vect();
        // Perp is the transverse difference in radial coordinate system (so
        // defined positive)
        Double_t ppv1_2 = p1.Perp(pg);
        hPperpGformtime->Fill(Gformtime, ppv1_2);
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
  hPtcharmvsGluonE->Write();
  TProfile *pPtcharmvsGluonE =
      (TProfile *)hPtcharmvsGluonE->ProfileX("pPtcharmvsGluonE");
  pPtcharmvsGluonE->GetYaxis()->SetTitle("p_{T} charm  (GeV)");
  pPtcharmvsGluonE->GetXaxis()->SetTitle("E_{g} (GeV)");
  pPtcharmvsGluonE->Write();
  hPt1Pt2->Write();
  TProfile *pPt1Pt2 = (TProfile *)hPt1Pt2->ProfileX("pPt1Pt2");
  pPt1Pt2->GetYaxis()->SetTitle("p_{T} HQ-1 (GeV)");
  pPt1Pt2->GetXaxis()->SetTitle("p_{T} HQ-2 (GeV)");
  pPt1Pt2->Write();
  hInvmassvsGluonE->Write();
  hInvmass->Write();
  TProfile *pGformtimevsGluonE =
      (TProfile *)hGformtimevsGluonE->ProfileX("pGformtimevsGluonE");
  hGformtimevsGluonE->Write();
  pGformtimevsGluonE->GetYaxis()->SetTitle("Formation time (fm/c)");
  pGformtimevsGluonE->GetXaxis()->SetTitle("E_{g} (GeV)");
  pGformtimevsGluonE->Write();
  hPperpGformtime->Write();
  TProfile *pPperpGformtime =
      (TProfile *)hPperpGformtime->ProfileX("pPperpGformtime");
  pPperpGformtime->GetYaxis()->SetTitle("c-gluon <q_{T}>");
  pPperpGformtime->GetXaxis()->SetTitle("Formation time (fm/c)");
  pPperpGformtime->Write();
  hPtcharmvsGformtime->Write();
  TProfile *pPtcharmvsGformtime =
      (TProfile *)hPtcharmvsGformtime->ProfileX("pPtcharmvsGformtime");
  pPtcharmvsGformtime->Write();
  hPxConservat->Write();
  fout->Write();
  return 0;
}
