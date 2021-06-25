// This is a simple macro to generate pythia predictions for the study 
// of gluon splitting processes. The basic idea is to find g->ccbar 
// topologies, identify the kinematic of the gluon and of the c/cbar
// and estimate the formation time of the gluon. In the future, a 
// simple parametrization of the broadening process will be added.

#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TList.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "fastjet/PseudoJet.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/ClusterSequenceArea.hh"
#include <ctime>
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <valarray>
#include <time.h>       /* time */
//#include <yaml.h>
//include <stdio.h>
//include <glib.h>
#include <yaml-cpp/yaml.h>

using namespace Pythia8;

int main(int argc, char* argv[]) {
    (void)argc;
    std::string mycase = argv[1];
    int cislo = -1;                 //unique number for each file
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

    //END OF CONFIGURATION

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

    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d",cislo));
    if (extramode=="mode2") {
        std::cout<<"Running with mode2"<<std::endl;
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
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation =5");
    }
    pythia.init();

    TFile *fout = new TFile(outputfile.data(), "recreate");

    fout->cd();

    // deltaR between the c-cbar pair defined in terms of eta and phi
    TH1F*hptdiff = new TH1F("hptdiff", ";#Delta pt (c#bar{c});Entries", 1000., 0, 3.);
    TH1F*hdeltaR = new TH1F("hdeltaR", ";#Delta R (c#bar{c});Entries", 1000., 0, 3.);
    TH1F*hetaG = new TH1F("hetaG", ";#eta(g);Entries", 1000., -5, 5.);
    // deltaR between the c-cbar pair defined in terms of eta and phi as a function of the gluon energy
    TH2F*hdeltaRvsGluonE = new TH2F("hdeltaRvsGluonE", ";E(g);#Delta R (c#bar{c})", 1000, 0., 100., 100., 0, 3.);
    // scatter plot pt charm vs pt anticharm
    TH2F*hPt1Pt2 = new TH2F("hPt1Pt2", ";p_{T}(c);p_{T}(cbar)", 200, 0., 200., 200., 0, 200.);
    TH2F*hPtcharmvsGluonE = new TH2F("hPtcharmvsGluonE", ";E(g);p_{T}(c)", 200, 0., 200., 200., 0, 200.);
    TH1F*hdeltaPhi = new TH1F("hdeltaPhi", ";#Delta #phi (c#bar{c});Entries", 1000., -5., 5.);
    TH2F*hQqbarmassvsGluonE = new TH2F("hQqbarmassvsGluonE", ";E(g);mass(c#bar{c})", 200, 0., 200., 200., 0, 200.);
    TH1F*hQqbarmass = new TH1F("hQqbarmass", ";mass (c#bar{c});Entries", 100., 0, 20.);
    TH2F*hQqbarformtimevsGluonE = new TH2F("hQqbarformtimevsGluonE", ";gluon energy; formation time(fm/c)", 200, 0., 200., 200., 0, 200.);
    TH2F*hqtvsformtime = new TH2F("hqtvsformtime", ";formation time (fm/c); q_{perp} (GeV)", 100, 0., 10., 100., 0, 20.);
    TH2F*hptcharmvsformtime = new TH2F("hptcharmvsformtime", ";formation time (fm/c); p_{T} charm (GeV)", 10, 0., 10., 100., 0, 20.);

    // Begin event loop. Generate event. Skip if error. List first one.
    for (int iEvent = 0; iEvent < maxnevents; ++iEvent) {
        pythia.next();
        for (int i = 0; i < pythia.event.size(); ++i) {
            double eCalc = pythia.event[i].eCalc();
            if (abs(eCalc/pythia.event[i].e() - 1.) > 1e-3) cout << " e mismatch, i = "
            << i << " e_nominal = " << pythia.event[i].e() << " e-from-p = "
            << eCalc << " m-from-e " << pythia.event[i].mCalc() << "\n";
            if(pythia.event[i].pT()<0 || pythia.event[i].pT()>1.e+5) continue;
            if(pythia.event[i].id()==hqpdg) {
		int mothercharmindex = pythia.event[i].mother1();
                int pdgmother = pythia.event[mothercharmindex].id();
   	        if (pdgmother!=21) continue;
		int daughter1 = pythia.event.daughterList(mothercharmindex)[0];
                int daughter2 = pythia.event.daughterList(mothercharmindex)[1];
	        int pdgdaughter1 = pythia.event[daughter1].id();		
	        int pdgdaughter2 = pythia.event[daughter2].id();
	        if (!((pdgdaughter1 == -hqpdg && pdgdaughter2 == hqpdg) || (pdgdaughter1 == hqpdg && pdgdaughter2 == -hqpdg))) continue; 
		int ndaughters = daughter2 - daughter1 + 1;
		if (ndaughters!=2) continue;
		if (pythia.event[mothercharmindex].e()<minEgluon || std::abs(pythia.event[mothercharmindex].eta())>1.) continue;
		double pt1 = pythia.event[daughter1].pT();
		double pt2 = pythia.event[daughter2].pT();
		double eta1 = pythia.event[daughter1].eta();
		double eta2 = pythia.event[daughter2].eta();
		double phi1 = pythia.event[daughter1].phi();
		double phi2 = pythia.event[daughter2].phi();
                double r = sqrt((eta2-eta1)*(eta2-eta1) + (phi2-phi1)*(phi2-phi1));
		hptdiff->Fill(pythia.event[daughter1].px() + pythia.event[daughter2].px() - pythia.event[mothercharmindex].px());
		hetaG->Fill(pythia.event[mothercharmindex].eta());
		hdeltaR->Fill(r);
		hdeltaPhi->Fill(phi2-phi1);
		hdeltaRvsGluonE->Fill(pythia.event[mothercharmindex].e(), r);
		hPtcharmvsGluonE->Fill(pythia.event[mothercharmindex].e(), pt1);
		hPt1Pt2->Fill(pt1, pt2);
		double massccbar = (pythia.event[daughter1].p() + pythia.event[daughter2].p()).mCalc();
                hQqbarmassvsGluonE->Fill(pythia.event[mothercharmindex].e(), massccbar);
		hQqbarmass->Fill(massccbar);
                // hcut*c=197.3MeV fm
		// In NU, 1 = 0.2 GeV * fm
		// 5 GeV-1 = fm 
		// GeV-1 = 0.2 fm
		// X GeV-1 = 0.2 X fm
		double formtime = 0.2 * pythia.event[daughter1].eCalc()/(2*pythia.event[daughter1].m2Calc());
                hQqbarformtimevsGluonE->Fill(pythia.event[mothercharmindex].eCalc(), formtime);
		TLorentzVector vg(pythia.event[mothercharmindex].eCalc(), pythia.event[mothercharmindex].px(), pythia.event[mothercharmindex].py(), pythia.event[mothercharmindex].pz());
		TLorentzVector v1(pythia.event[daughter1].e(), pythia.event[daughter1].px(), pythia.event[daughter1].py(), pythia.event[daughter1].pz());
		TLorentzVector v2(pythia.event[daughter2].e(), pythia.event[daughter2].px(), pythia.event[daughter2].py(), pythia.event[daughter2].pz());
		TVector3 pg = vg.Vect();
		TVector3 p1 = v1.Vect();
		TVector3 p2 = v2.Vect();
		Double_t ppv1_2 = p1.Perp(pg);
		hqtvsformtime->Fill(formtime, ppv1_2);
		hptcharmvsformtime->Fill(formtime, pt1);
	    }
            // End of event loop. Statistics. Histogram. Done.
        }
    }
    pythia.stat();
    printf("nAccepted %ld, nTried %ld\n", pythia.info.nAccepted(), pythia.info.nTried());
    printf("pythia.info.sigmaGen() %f\n", pythia.info.sigmaGen());
    hetaG->Write();
    hdeltaR->Write();
    hdeltaPhi->Write();
    hdeltaRvsGluonE->Write();
    hPtcharmvsGluonE->Write();
    TProfile * pPtcharmvsGluonE = (TProfile*)hPtcharmvsGluonE->ProfileX("pPtcharmvsGluonE");
    pPtcharmvsGluonE->Write();
    hPt1Pt2->Write();
    hQqbarmassvsGluonE->Write();
    hQqbarmass->Write();
    TProfile * pQqbarformtimevsGluonEg = (TProfile*)hQqbarformtimevsGluonE->ProfileX("pQqbarformtimevsGluonE");
    hQqbarformtimevsGluonE->Write();
    pQqbarformtimevsGluonEg->GetYaxis()->SetTitle("Formation time (fm/c)"); 
    pQqbarformtimevsGluonEg->GetXaxis()->SetTitle("E_{g} (GeV)"); 
    pQqbarformtimevsGluonEg->Write();
    hqtvsformtime->Write();
    TProfile * pqtvsformtime = (TProfile*)hqtvsformtime->ProfileX("pqtvsformtime");
    pqtvsformtime->GetYaxis()->SetTitle("c-gluon <q_{T}>"); 
    pqtvsformtime->GetXaxis()->SetTitle("Formation time (fm/c)"); 
    pqtvsformtime->Write();
    hptcharmvsformtime->Write();
    TProfile * pptcharmvsformtime = (TProfile*)hptcharmvsformtime->ProfileX("pptcharmvsformtime");
    pptcharmvsformtime->Write();
    hptdiff->Write();
    fout->Write();
    return 0;
}
