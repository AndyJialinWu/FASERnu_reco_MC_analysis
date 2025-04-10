#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCut.h"

#include "../DiscRecoMC/include/Fedra2022AnaRecoMC.h"

void plot_data(){
    
    gStyle->SetStatX(0.85);
    gStyle->SetTitleOffset(0.9,"X");
    
    gStyle->SetTitleSize(0.05,"XYZ");
    gStyle->SetLabelSize(0.04,"XYZ");
    
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);

    std::vector<std::string> *EvtNClist[3];
    
    for(int mcIt=0; mcIt<3; mcIt++){
        EvtNClist[mcIt] = new std::vector<std::string>();  // Allocate each vector
        ReadEventList(EvtNClist[mcIt], EventListFile[mcIt]);
        std::cout<<"# NC events from "<<EventListFile[mcIt]<<" = "<<EvtNClist[mcIt]->size()<<std::endl;
    }
    //PrintEventList(EvtNClist[0]);

    // MC 200025 200026 200035 Loop
    for(int mcIt=0; mcIt<3; mcIt++){
    
        int NEvts = EvtNClist[mcIt]->size();
        //double EventWeight = Run3ExpNC[mcIt] / NEvts;

        // Event Loop
        for(int EvtIt=0; EvtIt<NEvts; EvtIt++){

            TString LTAddress = EventDataDir[mcIt] + EvtNClist[mcIt]->at(EvtIt);
            TString LTname = LTAddress + "/linked_tracks.root";
            TFile *f_trfile = new TFile(LTname, "READ");
            TTree *tracks = (TTree *)f_trfile->Get("tracks");
            double xmax = tracks->GetMaximum("s.eX");
            double xmin = tracks->GetMinimum("s.eX");
            double ymax = tracks->GetMaximum("s.eY");
            double ymin = tracks->GetMinimum("s.eY");
            double pidmax = tracks->GetMaximum("s.ePID");
            double pidmin = tracks->GetMinimum("s.ePID");
            
            printf("range = %.0f %.0f %.0f %.0f  pid %.0f - %.0f\n", xmin, xmax, ymin, ymax, pidmin, pidmax);
            
            double xc = (xmax+xmin)/2.0;
            double yc = (ymax+ymin)/2.0;
            
            double range = 4000;
            
            int npl = (int) pidmax - (int) pidmin + 1;
            
            TCut cut = Form("nseg>=3&&sqrt(t.eTX**2+t.eTY**2)<0.1&&abs(s.eX-%.0f)<%.0f&&abs(s.eY-%.0f)<%.0f", xc, range, yc, range);
            TCanvas *c1 = new TCanvas("c1");
            c1->Divide(2,2);
            int ic1 = 1;
            
            c1->cd(ic1++);
            
            TH2F *hxy = new TH2F("hxy","",200, xc-range, xc+range, 200, yc-range, yc+range);
            tracks->Draw("s.eY:s.eX >>hxy",cut,"colz");
            double wx = hxy->GetXaxis()->GetBinWidth(1);
            double wy = hxy->GetYaxis()->GetBinWidth(1);
            hxy->SetXTitle("X");
            hxy->SetYTitle("Y");
            hxy->SetZTitle("segments/cm^{2}/film");
            
            hxy->Scale(1e8/wx/wy/npl);
            hxy->Draw("colz");
            
            c1->cd(ic1++);
            TH2F *htxty = new TH2F("htxty","",1000,-0.1,0.1, 1000,-0.1,0.1);
            tracks->Draw("s.eTY:s.eTX >>htxty",cut,"colz");
            htxty->SetXTitle("tan#thetax");
            htxty->SetYTitle("tan#thetay");
            htxty->SetZTitle("segments/cm^{2}/mrad^{2}/film");
            
            double wtx = htxty->GetXaxis()->GetBinWidth(1);
            double wty = htxty->GetYaxis()->GetBinWidth(1);
            
            htxty->Scale(1e-6/wtx/wty/npl*1e8/range/range/2/2);
            htxty->Draw("colz");
            
            c1->cd(ic1++);
            
            Int_t ix,iy, iz;
            htxty->GetMaximumBin(ix, iy, iz);
            TH2F *htxty2 = (TH2F *) htxty->Clone("htxty2");
            htxty2->GetXaxis()->SetRange(ix-100,ix+100);
            htxty2->GetYaxis()->SetRange(iy-100,iy+100);
            htxty2->Draw("colz");
            
            c1->cd(ic1++);
            TH2F *htxty3 = (TH2F *) htxty2->Clone("htxty3");
            htxty3->SetMaximum(250);
            htxty3->Draw("colz");
            
            c1->Print(LTAddress + "/segment_density.pdf");
            
            TFile f(LTAddress + "/segment_density.root","recreate");
            ((TH2F *)hxy->Clone("hXY"))->Write();
            ((TH2F *)htxty->Clone("hTXTY"))->Write();
            f.Close();

            delete c1;
            delete hxy; 
            delete htxty;
            delete htxty2;
            delete htxty3;
            f_trfile->Close();
            delete f_trfile;

        }// Event Loop
    
    }// MC 200025 200026 200035 Loop    
    
}
