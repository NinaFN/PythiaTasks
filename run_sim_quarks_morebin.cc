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
    TFile* outFile = new TFile("sim_tuples_quarks_pairprod.root", "RECREATE");

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
    vector<TNtuple*> quarkTuples(nBins);
    

    vector<double> binLuminocity(nBins); // luminocity from generated process sigma to calculate cross-sections
    vector<double> binEvCount(nBins);
    vector<double> weightSums(nBins);
    vector<double> sigmaGens(nBins);

    for (int i = 0; i < nBins; ++i) {
        //Decay Map:
            // 0: quark pair does not decay to RL or SL
            // 1: quark pair decays to mu
            // 2: quark pair decays to e
            // 3: quark pair decays to both e and mu
            // 4: quark pair decays to both e and mu, each lepton comes from seperate quark
            // 5: quark pair decays to both e and mu, each lepton comes from seperate quark, opp sign
            // 6: spicy 5 where both leptons are in the correct eta range

        //quarkTuples[i] = new TNtuple("quarks", "quarks", "binTag:event:id:pt:eta:decayMap1:decayMap2");        
        quarkTuples[i] = new TNtuple("quarks", "quarks", "binTag:event:id:pt1:ptHat:eta1:eta2:phi1:phi2:decayMap");        

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
            

            std::vector<double> elecEtas_q1;
            std::vector<double> elecEtas_q2;
            std::vector<double> muonEtas_q1;
            std::vector<double> muonEtas_q2;

            std::vector<double> elecCharges_q1;
            std::vector<double> elecCharges_q2;
            std::vector<double> muonCharges_q1;
            std::vector<double> muonCharges_q2;
            

            double pTHat  = pythia.info.pTHat();

            bool decayElec_q1 = false;
            bool decayElec_q2 = false;
            bool decayMuon_q1 = false;
            bool decayMuon_q2 = false;

            bool barrelHit = false;
            bool forwardHit = false;

            if(pythia.info.isNonDiffractive()){code = pythia.info.codeSub();}
            else{code = pythia.info.code();}

            if (softQCD && iBin < 2 && pythia.info.isNonDiffractive()
            && pTHat > binEdges[iBin+1]) continue;

            if (pTHat < binEdges[iBin]) continue;

            hardPtPart->Fill(pTHat);
            eventCount++;
            
            //cout << "====START OF NEW EVENT====" << endl;
            
            if(121<=code && code<=124) //only checks quarks if they come from a hardest process which produces a heavy flavor pair
            {
                for (int i = 0; i < pythia.event.size(); ++i) {
                    //cout <<"---------------------------------------- New Event"<<endl;
                    decayMap = 0;
                    //decayMap2 = 0;
                    
                    int particleID = pythia.event[i].id();
                    int particleStatus = pythia.event[i].status();
                    double particlePt = pythia.event[i].pT();
                    double particleEta = pythia.event[i].eta();
                    int particleCharge = pythia.event[i].charge();
                    

                    if (particleStatus==-23)
                    {
                        /*cout<<endl;
                        cout<<"EVENT: "<<iEvent<<endl;
                        cout<<"hardest process code: "<<code<<endl;
                        cout<<"particle from hardest process: "<<particleID<<endl;
                        cout<<"next particle in event record: "<<pythia.event[i+1].id()<<endl;
                        cout<<"status: "<<particleStatus<<endl<<endl;

                        cout<<"Quark 1:"<<endl;*/
                        for(int child: pythia.event[i].daughterListRecursive())
                        {
                            if(std::abs(pythia.event[child].id())==11)
                            {   
                                decayElec_q1 = true;
                                elecEtas_q1.push_back(pythia.event[child].eta());
                                elecCharges_q1.push_back(pythia.event[child].charge());
        
                                //cout<<"child particle found: "<<pythia.event[child].id()<<endl;
                            }

                            if(std::abs(pythia.event[child].id())==13)
                            {   
                                decayMuon_q1 = true;
                                muonEtas_q1.push_back(pythia.event[child].eta());
                                muonCharges_q1.push_back(pythia.event[child].charge());

                                //cout<<"child particle found: "<<pythia.event[child].id()<<endl;
                            }
                        }
                        //cout<<"Quark 2:"<<endl;
                        for(int child: pythia.event[i+1].daughterListRecursive())
                        {
                            if(std::abs(pythia.event[child].id())==11)
                            {
                                decayElec_q2 = true;
                                elecEtas_q2.push_back(pythia.event[child].eta());
                                elecCharges_q2.push_back(pythia.event[child].charge());
                                
                                //cout<<"child particle found: "<<pythia.event[child].id()<<endl;
                            }

                            if(std::abs(pythia.event[child].id())==13)
                            {
                                decayMuon_q2 = true;
                                muonEtas_q2.push_back(pythia.event[child].eta());
                                muonCharges_q2.push_back(pythia.event[child].charge());
                                
                                //cout<<"child particle found: "<<pythia.event[child].id()<<endl;
                            }
                        }

                        if(decayMuon_q1==true || decayMuon_q2==true){decayMap = 1;}
                        if(decayElec_q1==true || decayElec_q2==true){decayMap = 2;}
                        if((decayElec_q1==true || decayElec_q2==true) && (decayMuon_q1==true || decayMuon_q2==true)){decayMap = 3;}
                        
                        if((decayElec_q1==true && decayMuon_q2==true))
                        {
                            decayMap=4;

                            for(double elecCharge: elecCharges_q1){
                                for (int muonCharge: muonCharges_q2){
                                    if(elecCharge*muonCharge == -1){
                                        decayMap=5;
                                        //cout<<"!! e-mu pair found?"<<endl;
                                        break;
                                    }
                                }
                            }

                            for(double elecEta: elecEtas_q1){if(elecEta <= 0.9 && elecEta>=-0.9){barrelHit=true;}}
                            for(double muonEta: muonEtas_q2){if(muonEta >= -4.0 && muonEta<=-2.5){forwardHit=true;}}
                            if(barrelHit==true && forwardHit==true){decayMap=6;}
                        }

                        else if(decayElec_q2==true && decayMuon_q1==true)
                        {
                            decayMap=4;

                            for(double elecCharge: elecCharges_q2){
                                for (int muonCharge: muonCharges_q1){
                                    if(elecCharge*muonCharge == -1){
                                        decayMap=5;
                                        //cout<<"!! e-mu pair found?"<<endl;
                                        break;
                                    }
                                }
                            }

                            for(double elecEta: elecEtas_q2){if(elecEta <= 0.9 && elecEta>=-0.9){barrelHit=true;}}
                            for(double muonEta: muonEtas_q1){if(muonEta >= -4.0 && muonEta<=-2.5){forwardHit=true;}}
                            if(barrelHit==true && forwardHit==true){decayMap=6;}
                        }


                        //pythia.event.list(true);

                        //binTag:event:id:pt1:pt2:eta1:eta2:phi1:phi2:decayMap
                        
                        quarkTuples[iBin]->Fill(iBin,iEvent,particleID,particlePt,pTHat,particleEta,pythia.event[i+1].eta(),pythia.event[i].phi(),pythia.event[i+1].phi(),decayMap);                    
                        break;
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
        quarkTuples[i]->Write(Form("quark%d", i));
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
