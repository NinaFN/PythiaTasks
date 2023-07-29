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
    TFile* outFile = new TFile("sim_tuples_leps_individual.root", "RECREATE");

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

        qmTuples[i] = new TNtuple("quark_muon", "quark_muon", "binTag:event:parentFlag:pairPartner:ptHat:idQ:idL:ptQ:ptL:etaQ:etaL:phiQ:phiL");        
        qeTuples[i] = new TNtuple("quark_elec", "quark_elec", "binTag:event:parentFlag:pairPartner:ptHat:idQ:idL:ptQ:ptL:etaQ:etaL:phiQ:phiL");        

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
            genEvents = N_events*30;
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

            std::vector<vector<double>> quarkElecPos;
            std::vector<vector<double>> quarkMuonPos;

            
            //cout << "====START OF NEW EVENT====" << endl;
            
            if(121<=code && code<=124) //only checks quarks if they come from a hardest process which produces a heavy flavor pair
            {
                for (int i = 0; i < pythia.event.size(); ++i) {
                    //cout <<"---------------------------------------- New Event"<<endl;
                    
                    int particleStatus = pythia.event[i].status();
                    int particleID = pythia.event[i].id();
                    double particlePt = pythia.event[i].pT();
                    double particleEta = pythia.event[i].eta();
                    double particlePhi = pythia.event[i].phi();
                    
                    if (particleStatus==-23)
                    {
                        //quarkCount++;
                        for(int child: pythia.event[i].daughterListRecursive())
                        {
                            if(std::abs(pythia.event[child].id())==11){quarkElecPos.push_back({(double)child, (double)particleID, particlePt, particleEta, particlePhi});}
                            if(std::abs(pythia.event[child].id())==13){quarkMuonPos.push_back({(double)child, (double)particleID, particlePt, particleEta, particlePhi});}
                        }
                        
                    }

                    if(particleStatus>0 && std::abs(particleID)==11)
                    {
                        int quarkFlagElec = -1;
                        int parentID = 0;
                        double parentPt = 0;
                        double parentEta = 0;
                        double parentPhi = 0;
                        int muonFlag = -1;

                        for(vector<double> child: quarkElecPos)
                        {
                            if((int)child[0]==i)
                            {
                                quarkFlagElec = 1;
                                parentID = (int)child[1];
                                parentPt = (int)child[2];
                                parentEta = child[3];
                                parentPhi = child[4];

                                if (!quarkMuonPos.empty()) {muonFlag=1;}
                                //cout<<"!! child elec found"<<endl;
                                break;
                            }
                        }

                        qeTuples[iBin]->Fill(iBin,iEvent,quarkFlagElec,muonFlag,pTHat,parentID,particleID,parentPt,particlePt,parentEta,particleEta,parentPhi,particlePhi);
                    }

                    if(particleStatus>0 && std::abs(particleID)==13)
                    {
                        int quarkFlagMuon = -1;
                        int parentID = 0;
                        double parentPt = 0;
                        double parentEta = 0;
                        double parentPhi = 0;
                        int elecFlag = -1;

                        for(vector<double> child: quarkMuonPos)
                        {
                            if((int)child[0]==i)
                            {
                                quarkFlagMuon = 1;
                                parentID = (int)child[1];
                                parentPt = (int)child[2];
                                parentEta = child[3];
                                parentPhi = child[4];

                                if (!quarkElecPos.empty()) {elecFlag=1;}
                                //cout<<"!! child muon found"<<endl;
                                break;
                            }
                        }

                        qmTuples[iBin]->Fill(iBin,iEvent,quarkFlagMuon,elecFlag,pTHat,parentID,particleID,parentPt,particlePt,parentEta,particleEta,parentPhi,particlePhi);
                    }
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
