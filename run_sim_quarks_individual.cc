#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;


int main() {
    // Turn SoftQCD on/off
    bool softQCD = true;

    // Pythia object
    Pythia pythia;

    // ROOT file for histograms
    TFile* outFile = new TFile("sim_tuples_quarks_individual.root", "RECREATE");

    // pTHat bins
    int nBins;
    const double* binEdges;
    const int* binWidths;
    if (softQCD) {
        nBins = 8;
        static const double tempArray[9] = {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0};
        //static const int tempArrayW[6] = {16,10,10,14,20,30};

        binEdges = &tempArray[0];
        //binWidths = &tempArrayW[0];
    } else {
        nBins = 5;
        static const double tempArray[6] = {16.0, 26.0, 36.0, 50.0, 70.0, 100.0};
        //static const int tempArrayW[5] = {10,10,14,20,30};

        binEdges = &tempArray[0];
        //binWidths = &tempArrayW[0];
    }


    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 120, 0.0, 120.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 120, 0.0, 120.0);

    // HF Cross Sections
    vector<TNtuple*> qeTuples(nBins);
    vector<TNtuple*> qmTuples(nBins);
    

    vector<double> binLuminocity(nBins); // luminocity from generated process sigma to calculate cross-sections
    vector<double> binEvCount(nBins);
    vector<double> weightSums(nBins);
    vector<double> sigmaGens(nBins);

    for (int i = 0; i < nBins; ++i) {

        qmTuples[i] = new TNtuple("quark_muon", "quark_muon", "binTag:event:decayFlag:etaDecFlag:ptHat:idQ:idL:ptQ:ptL:etaQ:etaL:phiQ:phiL");        
        qeTuples[i] = new TNtuple("quark_elec", "quark_elec", "binTag:event:decayFlag:etaDecFlag:ptHat:idQ:idL:ptQ:ptL:etaQ:etaL:phiQ:phiL");        

    }

    // Number of events to generate per bin.
    int N_events = 50000;
    

    int decayMap = 0;
    //int decayMap2 = 0;

    int code = 0;


    int genEvents;

    // run events for each ptHat bin 
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (softQCD && iBin < 1) {
            pythia.readString("HardQCD:all = off");
            pythia.readString("HardQCD:hardccbar = off");
            pythia.readString("HardQCD:hardbbbar = off");
            pythia.readString("SoftQCD:nonDiffractive = on");
            genEvents = N_events*10;
        } else {
            // set pythia initialization variables
            genEvents = N_events;
            pythia.readString("HardQCD:all = off");
            pythia.readString("HardQCD:hardccbar = on");
            pythia.readString("HardQCD:hardbbbar = on");
            pythia.readString("SoftQCD:nonDiffractive = off");
        }

        pythia.readString("Beams:eCM = 13700.");
        pythia.readString("Tune:pp = 14");
        // pythia.readString("411:onMode=off");
        // pythia.readString("411:onIfAny=13");
        pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
        pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);
        pythia.init();

        int eventCount = 0;

        hardPtPart->Reset();
        cout<<"--------------------New Bin--------------------"<<endl;

        for (int iEvent = 0; iEvent < genEvents; ++iEvent) {
            //cout<<"\nEvent No: "<<iEvent+1<<endl;

            if (!pythia.next()) {continue;}
            
            

            double pTHat  = pythia.info.pTHat();

            if(pythia.info.isNonDiffractive()){code = pythia.info.codeSub();}
            else{code = pythia.info.code();}

            if (softQCD && iBin < 1 && pythia.info.isNonDiffractive()
            && pTHat > binEdges[iBin+1]) continue;

            if (pTHat < binEdges[iBin]) continue;

            hardPtPart->Fill(pTHat);
            eventCount++;

            int decayElec = -1;
            int decayMuon = -1;
            int barrelHit = -1;
            int forwardHit = -1;

            int elecID = 0;
            double elecPt = 0;
            double elecEta = 0;
            double elecPhi = 0;

            int muonID = 0;
            double muonPt = 0;
            double muonEta = 0;
            double muonPhi = 0;

            int quarkCount = 0;

            
            //cout << "====START OF NEW EVENT====" << endl;
            
            if(121<=code && code<=124) //only checks quarks if they come from a hardest process which produces a heavy flavor pair
            {
                for (int i = 0; i < pythia.event.size(); ++i) {
                    //cout <<"---------------------------------------- New Event"<<endl;
                    
                    int particleStatus = pythia.event[i].status();
                    
                    if (particleStatus==-23)
                    {
                        quarkCount++;

                        int particleID = pythia.event[i].id();
                        double particlePt = pythia.event[i].pT();
                        double particleEta = pythia.event[i].eta();
                        double particlePhi = pythia.event[i].phi();

                        for(int child: pythia.event[i].daughterListRecursive())
                        {
                            if(std::abs(pythia.event[child].id())==11)
                            {
                                decayElec = 1;
                            
                                elecID = pythia.event[child].id();
                                elecPt = pythia.event[child].pT();
                                elecEta = pythia.event[child].eta();
                                elecPhi = pythia.event[child].phi();

                                if(elecEta<=0.9 && elecEta>=-0.9)
                                {
                                    barrelHit=1;
                                    break; //prioritize storing info from leptons that do decay to correct eta range
                                }
                            }
                        }

                        for(int child: pythia.event[i].daughterListRecursive())
                        {
                            if(std::abs(pythia.event[child].id())==13)
                            {
                                decayMuon = 1;
                            
                                muonID = pythia.event[child].id();
                                muonPt = pythia.event[child].pT();
                                muonEta = pythia.event[child].eta();
                                muonPhi = pythia.event[child].phi();

                                if(muonEta>=-4.5 && muonEta<=-2.5)
                                {
                                    forwardHit=1;
                                    break; //prioritize storing info from leptons that do decay to correct eta range
                                }
                            }
                            
                        }
                        
                        //pythia.event.list(true);
                        
                        //binTag:event:decayFlag:etaDecFlag:ptHat:idQ:idL:ptQ:ptL:etaQ:etaL:phiQ:phiL
                        qeTuples[iBin]->Fill(iBin,iEvent,decayElec,barrelHit,pTHat,particleID,elecID,particlePt,elecPt,particleEta,elecEta,particlePhi,elecPhi);
                        qmTuples[iBin]->Fill(iBin,iEvent,decayMuon,forwardHit,pTHat,particleID,muonID,particlePt,muonPt,particleEta,muonEta,particlePhi,muonPhi);

                    }

                    if(quarkCount>2){break;}
                }
            }
        }

        // cross-section for the bin
        //double luminocity_hard = N_events/(pythia.info.sigmaGen()*pow(10,9));
        
        double luminocity_hard = (pythia.info.weightSum())/(pythia.info.sigmaGen()*pow(10,9));

        binLuminocity[iBin] = luminocity_hard;
        binEvCount[iBin] = eventCount;

        weightSums[iBin] = pythia.info.weightSum();
        sigmaGens[iBin] = pythia.info.sigmaGen();

        hardPtPart->Scale(1/luminocity_hard,"width");

        // add to final distribution
        hardPt->Add(hardPtPart);
    }

    for (int i = 0; i < nBins; ++i) 
    {   
        qeTuples[i]->Write(Form("qe%d", i));
        qmTuples[i]->Write(Form("qm%d", i));
    }

    outFile->WriteObject(&binLuminocity, "luminocities");
    outFile->WriteObject(&binEvCount, "eventCounts");
    outFile->WriteObject(&weightSums, "weightSums");
    outFile->WriteObject(&sigmaGens, "sigmaGens");
    
    // Total Cross Section
    TCanvas *canvasTotal = new TCanvas("total_sigma","total_sigma");
    gPad->SetLogy();

    hardPt->SetLineColor(1);
    hardPt->Draw();

    canvasTotal->Write();



    delete outFile;

    return 0;
}
