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
    TFile* outFile = new TFile("sim_tuples.root", "RECREATE");

    // pTHat bins
    int nBins;
    const double* binEdges;
    if (softQCD) {
        nBins = 7;
        static const double tempArray[8] = {0.0, 8.0, 16.0, 26.0, 36.0, 50.0, 70.0, 100.0};
        binEdges = &tempArray[0];
    } else {
        nBins = 5;
        static const double tempArray[6] = {16.0, 26.0, 36.0, 50.0, 70.0, 100.0};
        binEdges = &tempArray[0];
    }


    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 100, 0.0, 100.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 100, 0.0, 100.0);

    //Counts
    TH1F *charmCount = new TH1F("charm_count"," ;N_{particles}", 10, 0.0, 10);
    TH1F *bottomCount = new TH1F("bottom_count"," ;N_{particles}", 10, 0.0, 10);

    TH1F *elecCount = new TH1F("elec_count","Number of electrons produced per single collision;N_{particles}", 30, 0.0, 30);
    TH1F *elecCountBoth = new TH1F("elec_count_both","", 30, 0.0, 30);
    TH1F *elecCountHeavy = new TH1F("elec_count_heavy","", 30, 0.0, 30);
    TH1F *elecCountHeavyBoth = new TH1F("elec_count_heavy_both","", 30, 0.0, 30);

    TH1F *muonCount = new TH1F("muon_count","Number of muons produced per single collision;N_{particles}", 10, 0.0, 10);
    TH1F *muonCountBoth = new TH1F("muon_count_both","", 10, 0.0, 10);
    TH1F *muonCountHeavy = new TH1F("muon_count_heavy","", 10, 0.0, 10);
    TH1F *muonCountHeavyBoth = new TH1F("muon_count_heavy_both","", 10, 0.0, 10);

    // HF Cross Sections
    vector<TNtuple*> muonTuples(nBins);
    vector<TNtuple*> elecTuples(nBins);
    vector<TNtuple*> weightTuples(nBins);

    vector<double> binLuminocity(nBins); // luminocity from generated process sigma to calculate cross-sections


    //vector<TNtuple*> muonTuplesFull(nBins);
    //vector<TNtuple*> elecTuplesFull(nBins);

    for (int i = 0; i < nBins; ++i) {
        muonTuples[i] = new TNtuple("muon", "muon", "binTag:event:pt:eta:hadron:sign:flagBoth:flagOpp:flagHeavy");
        elecTuples[i] = new TNtuple("electron", "electron", "binTag:event:pt:eta:hadron:sign:flagBoth:flagOpp:flagHeavy");
        weightTuples[i] = new TNtuple("weight","weight","binTag:weight");
        
        //muonTuplesFull[i] = new TNtuple("muon", "muon", "pt:eta:hadron:hasBoth");
        //elecTuplesFull[i] = new TNtuple("electron", "electron", "pt:eta:hadron:hasBoth");
    }

    // Number of events to generate per bin.
    int N_events = 10000;

    bool hasElecPos;
    bool hasElecNeg;
    bool hasMuPos;
    bool hasMuNeg;

    int bothFlag;
    int heavyFlag;
    int oppFlag;

    int particleCount;



    // run events for each ptHat bin 
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (softQCD && iBin < 2) {
            pythia.readString("HardQCD:all = off");
            pythia.readString("SoftQCD:nonDiffractive = on");
        } else {
            // set pythia initialization variables
            pythia.readString("HardQCD:all = on");
            // pythia.readString("HardQCD:hardccbar = on");
            // pythia.readString("HardQCD:hardbbbar = on");
            pythia.readString("SoftQCD:nonDiffractive = off");
        }

        pythia.readString("Beams:eCM = 13700.");
        pythia.readString("Tune:pp = 14");
        // pythia.readString("411:onMode=off");
        // pythia.readString("411:onIfAny=13");
        pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
        pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);
        pythia.init();

        hardPtPart->Reset();
        cout<<"--------------------New Bin--------------------"<<endl;

        for (int iEvent = 0; iEvent < N_events; ++iEvent) {
            //cout<<"\nEvent No: "<<iEvent+1<<endl;

            if (!pythia.next()) {continue;}
            
            hasElecPos = false;
            hasElecNeg = false;
            hasMuPos = false;
            hasMuNeg = false;

            bothFlag = -1;
            heavyFlag = -1;
            oppFlag = -1;
            
            particleCount = 0;

            for (int i = 0; i < pythia.event.size(); ++i) 
            {

                if(pythia.event[i].status() > 0)
                {
                    if(pythia.event[i].id() == 13){hasMuPos = true;}
                    if(pythia.event[i].id() == -13){hasMuNeg = true;}
                    if(pythia.event[i].id() == 11){hasElecPos = true;}
                    if(pythia.event[i].id() == -11){hasElecNeg = true;}
                    
                }
            }


            if((hasMuPos == true || hasMuNeg == true)&&(hasElecPos == true || hasElecNeg == true))
            {
                bothFlag = 1;
                if((hasMuPos == true && hasElecNeg == true)||(hasMuNeg == true && hasElecPos == true))
                {
                    oppFlag = 1;
                }
            }

            double pTHat  = pythia.info.pTHat();

            if (softQCD && iBin < 2 && pythia.info.isNonDiffractive()
            && pTHat > binEdges[iBin+1]) continue;

            if (pTHat < binEdges[iBin]) continue;

            hardPtPart->Fill(pTHat);
            

            //cout << "====START OF NEW EVENT====" << endl;

            for (int i = 0; i < pythia.event.size(); ++i) {
                //cout <<"---------------------------------------- New Event";
                particleCount++;
                
                int particleID = abs(pythia.event[i].id());
                int particleStatus = pythia.event[i].status();
                int particlePt = pythia.event[i].pT();
                int particleEta = pythia.event[i].eta();
                int particleCharge = pythia.event[i].charge();
                //bool heavyConf = false;
                int heavyConf = -1;

                if (particleStatus > 0 && std::abs(particleID) == 13) 
                { // muon
                    int motherIndex = pythia.event[i].mother1();
                    int firstMotherID = abs(pythia.event[motherIndex].id());

                    //cout << "mu";
                    bool motherHadronStatus = pythia.event[motherIndex].isHadron();
                    int prevIndex = -1;
                    while (motherHadronStatus) {
                        //cout << " -> " << pythia.event[motherIndex].name();
                        prevIndex = motherIndex;
                        motherIndex = pythia.event[motherIndex].mother1();
                        motherHadronStatus = pythia.event[motherIndex].isHadron();
                    }

                    if (prevIndex != -1) {
                        int lastMotherID = abs(pythia.event[prevIndex].id());
                        //cout << " -> [";
                        for (int index: pythia.event[prevIndex].motherList()) {
                            if (!pythia.event[index].isGluon()) {
                                //cout << pythia.event[index].name() << " ";
                                cout<<" ";
                            }
                        }
                        //cout << "]" << endl;
                    } else {
                        int lastMotherID = firstMotherID;
                        //cout << " -> " << pythia.event[motherIndex].name() << endl;
                    }

                    muonTuples[iBin]->Fill(iBin,iEvent,particlePt, particleEta, firstMotherID, particleCharge, bothFlag, oppFlag, heavyFlag);                    
                }


                if (particleStatus > 0 && std::abs(particleID) == 11) 
                { // electron
                    int motherIndex = pythia.event[i].mother1();
                    int firstMotherID = abs(pythia.event[motherIndex].id());

                    //cout << "e";
                    bool motherHadronStatus = pythia.event[motherIndex].isHadron();
                    int prevIndex = -1;
                    while (motherHadronStatus) {
                        //cout << " -> " << pythia.event[motherIndex].name();
                        prevIndex = motherIndex;
                        motherIndex = pythia.event[motherIndex].mother1();
                        motherHadronStatus = pythia.event[motherIndex].isHadron();
                    }

                    if (prevIndex != -1) {
                        int lastMotherID = abs(pythia.event[prevIndex].id());
                        //cout << " -> [";
                        for (int index: pythia.event[prevIndex].motherList()) {
                            if (!pythia.event[index].isGluon()) {
                                cout<<" ";
                                //cout << pythia.event[index].name() << " ";
                            }
                        }
                        //cout << "]" << endl;
                    } else {
                        int lastMotherID = firstMotherID;
                        //cout << " -> " << pythia.event[motherIndex].name() << endl;
                    }

                    elecTuples[iBin]->Fill(iBin,iEvent, particlePt, particleEta, firstMotherID, particleCharge, bothFlag, oppFlag, heavyFlag);                    
                }
                //pythia.event.list(true);
            }
            //cout<<"\n Particle count: "<<particleCount<<endl;

        }

        // cross-section for the bin
        double luminocity_hard = N_events/(pythia.info.sigmaGen()*pow(10,9));
        binLuminocity[iBin] = luminocity_hard;

        weightTuples[iBin]->Fill(iBin,luminocity_hard);

        hardPtPart->Scale(1/luminocity_hard, "width");

        // add to final distribution
        hardPt->Add(hardPtPart);
    }

    for (int i = 0; i < nBins; ++i) 
    {   
        muonTuples[i]->Write(Form("muon%d", i));
        elecTuples[i]->Write(Form("elec%d", i));
    }

    outFile->WriteObject(&binLuminocity, "luminocities");
    //outFile->WriteObject(&muonTuples, "muons");
    //outFile->WriteObject(&elecTuples, "electrons");
    
    // Total Cross Section
    TCanvas *canvasTotal = new TCanvas("total_sigma","total_sigma");

    hardPt->SetLineColor(1);
    hardPt->Draw();

    canvasTotal->Write();



    delete outFile;

    return 0;
}