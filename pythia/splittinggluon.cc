// This is a simple macro to generate pythia predictions

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

    int maxnevents = node["maxneventsperjob"].as<int>();
    int tune = node["tune"].as<int>();
    int beamidA = node["beamidA"].as<int>();
    int beamidB = node["beamidB"].as<int>();
    float eCM = node["eCM"].as<float>();
    float ptmingluon = node["ptmingluon"].as<float>();
    float pTHatMin = node["pTHatMin"].as<float>();
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
    TH1F*hdeltaR = new TH1F("hdeltaR", ";#Delta R (c#bar{c});Entries", 1000., 0, 3.);
    TH2F*hdeltaRvsGluonPt = new TH2F("hdeltaRvsGluonPt", ";p_{T}(g);#Delta R (c#bar{c})", 50, 0., 50., 30., 0, 3.);
    TH1F*hdeltaPhi = new TH1F("hdeltaPhi", ";#Delta #phi (c#bar{c});Entries", 1000., -5., 5.);

    // Begin event loop. Generate event. Skip if error. List first one.
    for (int iEvent = 0; iEvent < maxnevents; ++iEvent) {
        pythia.next();
        for (int i = 0; i < pythia.event.size(); ++i) {
            if(pythia.event[i].pT()<0 || pythia.event[i].pT()>1.e+5) continue;
            if(pythia.event[i].id()==4) {
		int mothercharmindex = pythia.event[i].mother1();
                int pdgmother = pythia.event[mothercharmindex].id();
		if (pythia.event[mothercharmindex].pT() < ptmingluon) continue;
   	        if (pdgmother!=21) continue;
		int daughter1 = pythia.event.daughterList(mothercharmindex)[0];
                int daughter2 = pythia.event.daughterList(mothercharmindex)[1];
	        int pdgdaughter1 = pythia.event[daughter1].id();		
	        int pdgdaughter2 = pythia.event[daughter2].id();
	        if (!((pdgdaughter1 == -4 && pdgdaughter2 == 4) || (pdgdaughter1 == 4 && pdgdaughter2 == -4))) continue;           

		if (std::abs(pythia.event[daughter1].y()<1.)) continue;
                if (std::abs(pythia.event[daughter2].y()<1.)) continue;

		int ndaughters = daughter2 - daughter1 + 1;
		if (ndaughters!=2) continue;
		//std::cout<<"daughter1 index="<<daughter1<<std::endl;
		//std::cout<<"daughter2 index="<<daughter2<<std::endl;
		//std::cout<<"ndaughters="<<ndaughters<<std::endl;
		double eta1 = pythia.event[daughter1].eta();
		double eta2 = pythia.event[daughter2].eta();
		double phi1 = pythia.event[daughter1].phi();
		double phi2 = pythia.event[daughter2].phi();
                double r = sqrt((eta2-eta1)*(eta2-eta1) + (phi2-phi1)*(phi2-phi1));
		hdeltaR->Fill(r);
		hdeltaPhi->Fill(phi2-phi1);
		hdeltaRvsGluonPt->Fill(pythia.event[mothercharmindex].pT(), r);
	    }
            // End of event loop. Statistics. Histogram. Done.
        }
    }
    pythia.stat();
    printf("nAccepted %ld, nTried %ld\n", pythia.info.nAccepted(), pythia.info.nTried());
    printf("pythia.info.sigmaGen() %f\n", pythia.info.sigmaGen());
    hdeltaR->Write();
    hdeltaPhi->Write();
    hdeltaRvsGluonPt->Write();
    fout->Write();
    return 0;
}
