#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

void read_sim() {
    TFile *f = TFile::Open("sim_tuples.root","READ");
    TFile *outf =  new TFile("sim_hists.root", "RECREATE");
    
    vector<double> *binLuminocity;
    vector<TNtuple*> muonTuples(7); //should be count of bins but 7 is easier
    vector<TNtuple*> elecTuples(7);

    f->GetObject("luminocities",binLuminocity);

    int binCount = 0;

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

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        muonTuples[binCount] = (TNtuple*)f->Get(Form("muon%d", binCount));
        elecTuples[binCount] = (TNtuple*)f->Get(Form("elec%d", binCount));

        ////Create Histograms

        //Muons
        muonPtPart->Reset();
        muonTuples[binCount]->Draw("pt>>muon_pt_part");
        muonPtPart->Scale(1/(*it), "width");
        muonPtTotal->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[binCount]->Draw("pt>>muon_pt_part", "flagBoth == 1");
        muonPtPart->Scale(1/(*it), "width");
        muonPtBoth->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[binCount]->Draw("pt>>muon_pt_part", "flagOpp == 1");
        muonPtPart->Scale(1/(*it), "width");
        muonPtOpp->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[binCount]->Draw("pt>>muon_pt_part", "hadron >= 411 && hadron <= 435 || hadron >= 511 && hadron <= 545 || hadron == 15 || hadron == 443 || hadron == 553");
        muonPtPart->Scale(1/(*it), "width");
        muonPtPrim->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[binCount]->Draw("pt>>muon_pt_part", "flagBoth == 1 && (hadron >= 411 && hadron <= 435 || hadron >= 511 && hadron <= 545 || hadron == 15 || hadron == 443 || hadron == 553)");
        muonPtPart->Scale(1/(*it), "width");
        muonPtBothPrim->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[binCount]->Draw("pt>>muon_pt_part", "flagOpp == 1 && (hadron >= 411 && hadron <= 435 || hadron >= 511 && hadron <= 545 || hadron == 15 || hadron == 443 || hadron == 553)");
        muonPtPart->Scale(1/(*it), "width");
        muonPtOppPrim->Add(muonPtPart);

        //Muon Geometry
        muonPtPart->Reset();
        muonTuples[binCount]->Draw("pt>>muon_pt_part", "-4.0<eta<2.5");
        muonPtPart->Scale(1/(*it), "width");
        muonPtGeo->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[binCount]->Draw("pt>>muon_pt_part", "-4.0<eta<2.5 && flagBoth==1");
        muonPtPart->Scale(1/(*it), "width");
        muonPtGeoBoth->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuples[binCount]->Draw("pt>>muon_pt_part", "-4.0<eta<2.5 && flagOpp==1");
        muonPtPart->Scale(1/(*it), "width");
        muonPtGeoOpp->Add(muonPtPart);
      
        binCount++;
    }


    ////Plotting
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

    delete outf;

}