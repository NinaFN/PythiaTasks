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
    TFile* outFile = new TFile("mymain09_2.root", "RECREATE");

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

    double binLuminocity[nBins]; // luminocity from generated process sigma to calculate cross-sections

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
    int N_events = 50000;

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

                    muonTuples[iBin]->Fill(iEvent,particlePt, particleEta, firstMotherID, particleCharge, bothFlag, oppFlag, heavyFlag);                    
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

                    elecTuples[iBin]->Fill(iEvent, particlePt, particleEta, firstMotherID, particleCharge, bothFlag, oppFlag, heavyFlag);                    
                }
                //pythia.event.list(true);
            }
            //cout<<"\n Particle count: "<<particleCount<<endl;

        }

        // cross-section for the bin
        double luminocity_hard = N_events/(pythia.info.sigmaGen()*pow(10,9));
        binLuminocity[iBin] = luminocity_hard;

        hardPtPart->Scale(1/luminocity_hard, "width");

        // add to final distribution
        hardPt->Add(hardPtPart);
    }

    
    
    
    
    // Generate histograms from bin tuples
    // muon
    TH1F *muonPtTotal = new TH1F("muon_full","Produced Muon Cross-Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *muonPtPart = new TH1F("muon_pt_part","", 35, 0.0, 70.0);
    TH1F *muonPtPrim = new TH1F("muon_pt_prim","", 35, 0.0, 70.0);
    TH1F *muonPtBoth = new TH1F("muon_pt_both","", 35, 0.0, 70.0);
    TH1F *muonPtBothPrim = new TH1F("muon_pt_both_prim","", 35, 0.0, 70.0);
    TH1F *muonPtOpp = new TH1F("muon_pt_opp","", 35, 0.0, 70.0);
    TH1F *muonPtOppPrim = new TH1F("muon_pt_opp_prim","", 35, 0.0, 70.0);

    TH1F *muonPtGeo = new TH1F("muon_pt_geo","", 35, 0.0, 70.0);
    TH1F *muonPtGeoBoth = new TH1F("muon_pt_geo_both","", 35, 0.0, 70.0);
    TH1F *muonPtGeoOpp = new TH1F("muon_pt_geo_opp","", 35, 0.0, 70.0);



    // electron
    TH1F *elecPtTotal = new TH1F("elec_full","Produced Electron Cross-Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *elecPtPart = new TH1F("elec_pt_part","", 35, 0.0, 70.0);
    TH1F *elecPtPrim = new TH1F("elec_pt_prim","", 35, 0.0, 70.0);
    TH1F *elecPtBoth = new TH1F("elec_pt_both","", 35, 0.0, 70.0);
    TH1F *elecPtBothPrim = new TH1F("elec_pt_both_prim","", 35, 0.0, 70.0);
    TH1F *elecPtOpp = new TH1F("elec_pt_opp","", 35, 0.0, 70.0);
    TH1F *elecPtOppPrim = new TH1F("elec_pt_opp_prim","", 35, 0.0, 70.0);
   
    TH1F *elecPtGeo = new TH1F("elec_pt_geo","", 35, 0.0, 70.0);
    TH1F *elecPtGeoBoth = new TH1F("elec_pt_geo_both","", 35, 0.0, 70.0);
    TH1F *elecPtGeoOpp = new TH1F("elec_pt_geo_opp","", 35, 0.0, 70.0);

    for (int i = 0; i < nBins; ++i) {
        //muon
        muonPtPart->Reset();
        muonTuples[i]->Draw("pt>>muon_pt_part");
        muonPtPart->Scale(1/binLuminocity[i], "width");
        muonPtTotal->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[i]->Draw("pt>>muon_pt_part", "flagBoth == 1");
        muonPtPart->Scale(1/binLuminocity[i], "width");
        muonPtBoth->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[i]->Draw("pt>>muon_pt_part", "flagOpp == 1");
        muonPtPart->Scale(1/binLuminocity[i], "width");
        muonPtOpp->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[i]->Draw("pt>>muon_pt_part", "hadron >= 411 && hadron <= 435 || hadron >= 511 && hadron <= 545 || hadron == 15 || hadron == 443 || hadron == 553");
        muonPtPart->Scale(1/binLuminocity[i], "width");
        muonPtPrim->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[i]->Draw("pt>>muon_pt_part", "flagBoth == 1 && (hadron >= 411 && hadron <= 435 || hadron >= 511 && hadron <= 545 || hadron == 15 || hadron == 443 || hadron == 553)");
        muonPtPart->Scale(1/binLuminocity[i], "width");
        muonPtBothPrim->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[i]->Draw("pt>>muon_pt_part", "flagOpp == 1 && (hadron >= 411 && hadron <= 435 || hadron >= 511 && hadron <= 545 || hadron == 15 || hadron == 443 || hadron == 553)");
        muonPtPart->Scale(1/binLuminocity[i], "width");
        muonPtOppPrim->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[i]->Draw("pt>>muon_pt_part", "-4.0<eta<-2.5");
        muonPtPart->Scale(1/binLuminocity[i], "width");
        muonPtGeo->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[i]->Draw("pt>>muon_pt_part", "-4.0<eta<-2.5 && flagBoth==1");
        muonPtPart->Scale(1/binLuminocity[i], "width");
        muonPtGeoBoth->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[i]->Draw("pt>>muon_pt_part", "-4.0<eta<-2.5 && flagOpp==1");
        muonPtPart->Scale(1/binLuminocity[i], "width");
        muonPtGeoOpp->Add(muonPtPart);

        //electron
        elecPtPart->Reset();
        elecTuples[i]->Draw("pt>>elec_pt_part");
        elecPtPart->Scale(1/binLuminocity[i], "width");
        elecPtTotal->Add(elecPtPart);

        elecPtPart->Reset();
        elecTuples[i]->Draw("pt>>elec_pt_part", "flagBoth == 1");
        elecPtPart->Scale(1/binLuminocity[i], "width");
        elecPtBoth->Add(elecPtPart);

        elecPtPart->Reset();
        elecTuples[i]->Draw("pt>>elec_pt_part", "flagOpp == 1");
        elecPtPart->Scale(1/binLuminocity[i], "width");
        elecPtOpp->Add(elecPtPart);

        elecPtPart->Reset();
        elecTuples[i]->Draw("pt>>elec_pt_part", "hadron >= 411 && hadron <= 435 || hadron >= 511 && hadron <= 545 || hadron == 15 || hadron == 443 || hadron == 553");
        elecPtPart->Scale(1/binLuminocity[i], "width");
        elecPtPrim->Add(elecPtPart);

        elecPtPart->Reset();
        elecTuples[i]->Draw("pt>>elec_pt_part", "flagBoth == 1 && (hadron >= 411 && hadron <= 435 || hadron >= 511 && hadron <= 545 || hadron == 15 || hadron == 443 || hadron == 553)");
        elecPtPart->Scale(1/binLuminocity[i], "width");
        elecPtBothPrim->Add(elecPtPart);

        elecPtPart->Reset();
        elecTuples[i]->Draw("pt>>elec_pt_part", "flagOpp == 1 && (hadron >= 411 && hadron <= 435 || hadron >= 511 && hadron <= 545 || hadron == 15 || hadron == 443 || hadron == 553)");
        elecPtPart->Scale(1/binLuminocity[i], "width");
        elecPtOppPrim->Add(elecPtPart);

        elecPtPart->Reset();
        elecTuples[i]->Draw("pt>>elec_pt_part", "-0.9<eta<0.9");
        elecPtPart->Scale(1/binLuminocity[i], "width");
        elecPtGeo->Add(elecPtPart);

        elecPtPart->Reset();
        elecTuples[i]->Draw("pt>>elec_pt_part", "-0.9<eta<0.9 && flagBoth==1");
        elecPtPart->Scale(1/binLuminocity[i], "width");
        elecPtGeoBoth->Add(elecPtPart);

        elecPtPart->Reset();
        elecTuples[i]->Draw("pt>>elec_pt_part", "-0.9<eta<0.9 && flagOpp==1");
        elecPtPart->Scale(1/binLuminocity[i], "width");
        elecPtGeoOpp->Add(elecPtPart);

    }

    //Plotting

    // Total Cross Section
    TCanvas *canvasTotal = new TCanvas("total_sigma","total_sigma");

    hardPt->SetLineColor(1);
    hardPt->Draw();

    canvasTotal->Write();

    // Output Muons
    TCanvas *canvasMuon = new TCanvas("Muon_sigma","Muon_sigma");

    muonPtTotal->SetLineColor(kRed);
   //muonPtTotal->SetMarkerStyle(20);
    //muonPtTotal->SetStats(0);
    muonPtTotal->Draw();

    muonPtPrim->SetLineColor(kOrange);
    //muonPtPrim->SetMarkerStyle(22);
    muonPtPrim->SetStats(0);
    muonPtPrim->Draw("SAME");

    muonPtBoth->SetLineColor(kGreen);
    //muonPtValid->SetMarkerStyle(21);
    muonPtBoth->SetStats(0);
    muonPtBoth->Draw("SAME");

    muonPtBothPrim->SetLineColor(kAzure+8);
    muonPtBothPrim->SetStats(0);
    muonPtBothPrim->Draw("SAME");

    muonPtOpp->SetLineColor(kViolet);
    //muonPtValid->SetMarkerStyle(21);
    muonPtOpp->SetStats(0);
    muonPtOpp->Draw("SAME");   

    muonPtOppPrim->SetLineColor(kViolet+2);
    muonPtOppPrim->SetStats(0);
    muonPtOppPrim->Draw("SAME");


    auto legendMuon = new TLegend();
    legendMuon->AddEntry(muonPtTotal,"Total #mu yield","l");
    legendMuon->AddEntry(muonPtPrim,"#mu yield from heavy quark primary decays","l");
    legendMuon->AddEntry(muonPtBoth,"#mu yield from collisions with any e production","l");
    legendMuon->AddEntry(muonPtBothPrim,"#mu yield from heavy quark primary decays in collisions with any e production","l");
    legendMuon->AddEntry(muonPtOpp,"#mu yield from collisions with opp-sign e production","l");
    legendMuon->AddEntry(muonPtOppPrim,"#mu yield from heavy quark primary decays in collisions with opp sign e production","l");
    legendMuon->Draw("SAME");

    canvasMuon->Write();

    //Geo Muon

    TCanvas *canvasMuonGeo = new TCanvas("Muon_sigma_geo","Muon_sigma_geo");

    muonPtTotal->SetLineColor(kBlack);
   //muonPtTotal->SetMarkerStyle(20);
    //muonPtTotal->SetStats(0);
    muonPtTotal->Draw();

    muonPtGeo->SetLineColor(kRed);
    muonPtGeo->SetStats(0);
    muonPtGeo->Draw("SAME");

    muonPtGeoBoth->SetLineColor(kGreen);
    muonPtGeoBoth->SetStats(0);
    muonPtGeoBoth->Draw("SAME");

    muonPtGeoOpp->SetLineColor(kBlue);
    muonPtGeoOpp->SetStats(0);
    muonPtGeoOpp->Draw("SAME");

    auto legendMuonGeo = new TLegend();
    legendMuonGeo->AddEntry(muonPtTotal,"Total #mu yield","l");
    legendMuonGeo->AddEntry(muonPtGeo,"#mu yield in forward arm #eta range","l");
    legendMuonGeo->AddEntry(muonPtGeoBoth,"#mu yield in forward arm #eta range in collisions with any e production","l");
    legendMuonGeo->AddEntry(muonPtGeoOpp,"#mu yield in forward arm #eta range in collisions with opposite signed e","l");
    legendMuonGeo->Draw("SAME");

    canvasMuonGeo->Write();


    // Output Electrons
    TCanvas *canvasElec = new TCanvas("Elec_sigma","Elec_sigma");

    elecPtTotal->SetLineColor(kRed);
    //elecPtTotal->SetMarkerStyle(20);
    //muonPtTotal->SetStats(0);
    elecPtTotal->Draw();

    elecPtPrim->SetLineColor(kOrange);
    //elecPtPrim->SetMarkerStyle(22);
    elecPtPrim->SetStats(0);
    elecPtPrim->Draw("SAME");

    elecPtBoth->SetLineColor(kGreen);
    //elecPtValid->SetMarkerStyle(21);
    elecPtBoth->SetStats(0);
    elecPtBoth->Draw("SAME");

    elecPtBothPrim->SetLineColor(kAzure+8);
    //elecPtValidPrim->SetMarkerStyle(23);
    elecPtBothPrim->SetStats(0);
    elecPtBothPrim->Draw("SAME");

    elecPtOpp->SetLineColor(kViolet);
    elecPtOpp->SetStats(0);
    elecPtOpp->Draw("SAME");

    elecPtOppPrim->SetLineColor(kViolet+2);
    elecPtOppPrim->SetStats(0);
    elecPtOppPrim->Draw("SAME");

    auto legendElec = new TLegend();
    legendElec->AddEntry(elecPtTotal,"Total e yield","l");
    legendElec->AddEntry(elecPtPrim,"e yield from heavy quark primary decays","l");
    legendElec->AddEntry(elecPtBoth,"e yield from collisions with #mu production","l");
    legendElec->AddEntry(elecPtBothPrim,"e yield from heavy quark primary decays in collisions with #mu production","l");
    legendElec->AddEntry(elecPtOpp,"e yield from collisions with opp-sign #mu production","l");
    legendElec->AddEntry(elecPtOppPrim,"e yield from heavy quark primary decays in collisions with opp sign #mu production","l");
    legendElec->Draw("SAME");

    canvasElec->Write();
    
    //Geo Electron

    TCanvas *canvasElecGeo = new TCanvas("Elec_sigma_geo","Elec_sigma_geo");

    elecPtTotal->SetLineColor(kBlack);
    elecPtTotal->Draw();

    elecPtGeo->SetLineColor(kRed);
    elecPtGeo->SetStats(0);
    elecPtGeo->Draw("SAME");

    elecPtGeoBoth->SetLineColor(kGreen);
    elecPtGeoBoth->SetStats(0);
    elecPtGeoBoth->Draw("SAME");

    elecPtGeoOpp->SetLineColor(kBlue);
    elecPtGeoOpp->SetStats(0);
    elecPtGeoOpp->Draw("SAME");

    auto legendElecGeo = new TLegend();
    legendElecGeo->AddEntry(elecPtTotal,"Total e yield","l");
    legendElecGeo->AddEntry(elecPtGeo,"e yield in central barrel #eta range","l");
    legendElecGeo->AddEntry(elecPtGeoBoth,"e yield in central barrel #eta range in collisions with any #mu production","l");
    legendElecGeo->AddEntry(elecPtGeoOpp,"e yield in central barrel #eta range in collisions with opposite signed #mu","l");
    legendElecGeo->Draw("SAME");

    canvasElecGeo->Write();

    delete outFile;

    return 0;
}