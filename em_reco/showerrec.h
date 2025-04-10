#ifndef __CINT__
#include<EdbSegP.h>
#include<EdbEDA.h>
#endif


class ShowerData{
    public:
    EdbSegP *axis;
    double E0;
    int pid;

    double ShowerMaximum;
    double NsegSum;
    double NsegSumShowerMax;
    double NsegShowerMax;
    double Nbg;
    double NbgShowerMax;  

    double Erec;

    void Print(){
        axis->PrintNice();
        printf("E0 = %lf pid=%3d\n",E0, pid);
        printf("Shower maximum at %.0f with %.0f segments\n", ShowerMaximum, NsegShowerMax);
        printf("NsegSum = %.0f (%.0f seg around shower max)\n", NsegSum, NsegSumShowerMax);
        printf("Erec = %.1lf\n", Erec);
        
    }  
    void PrintData(){
        FILE *fp = fopen("shower.txt","wt");
        if(ferror(fp)) return;
        fprintf(fp,"%d %d %.1f %.1f %.1f %.4f %.4f %.1lf %.1lf %.1lf %.1lf %.0lf %.0lf\n",
            axis->MCEvt(), axis->MCTrack(), axis->X(), axis->Y(), axis->Z(), axis->TX(), axis->TY(),
            E0, Erec, ShowerMaximum, NsegShowerMax, NsegSum, NsegSumShowerMax);
        fclose(fp);
    }  
    void ReadData(char *filename="shower.txt"){
        FILE *fp = fopen(filename,"rt");
        char buf[512];
        if(ferror(fp)||feof(fp)) return;
        fgets(buf, sizeof(buf), fp);       
        fclose(fp);
        int iev;
        float x,y,z,tx,ty;
        sscanf(buf,"%d %d %f %f %f %f %f %lf %lf %lf %lf %lf\n",
            &iev, &pid, &x, &y, &z, &tx, &ty,
            &E0, &ShowerMaximum, &NsegShowerMax, &NsegSum, &NsegSumShowerMax);
        axis->Set(0,x,y,tx,ty,0,0);
        axis->SetZ(z);
        axis->SetMC(iev, pid);
        PrintData();
    } 
};


double ErecShowerMaximum(double showermax){
    // with d<100, dth<0.01, dmin<50
    // memo

    // nt->Draw("log10(E):ShowerMaximum","E>10&&index==0&&NsegSum>50&&E>50&&ShowerMaximum>=10&&ShowerMaximum<40","same");
    // TGraph *gr = (TGraph *) c1->FindObject("Graph");
    // if(gr) gr->Fit("pol1","ROB=0.9","",12,40);
    return TMath::Power(10,0.060*showermax+1.39);

}

double ErecNsegSumShowerMax(double nsegSumShowerMax){
      
    double p0                        =      22.9916;
    double p1                        =       2.4266;
    double p2                        =    0.0189621;
    double p3                        =  8.45026e-05;

    double x = nsegSumShowerMax;
    return p0 + p1*x + p2*TMath::Power(x,2) + p3*TMath::Power(x,3);

}

ShowerData shower;

void FitTrack(EdbTrackP *t){
    if(t->N()<=1) return;

    TGraph grx,gry;
    
    for(int i=0; i<t->N()&&i<=3; i++){
        EdbSegP *s = t->GetSegment(i);
        grx.SetPoint(grx.GetN(), s->Z(), s->X());
        gry.SetPoint(gry.GetN(), s->Z(), s->Y());
    }

    grx.Fit("pol1");
    gry.Fit("pol1");
    t->SetTX(grx.GetFunction("pol1")->GetParameter(1));
    t->SetTY(gry.GetFunction("pol1")->GetParameter(1));
    
    t->PrintNice();

}

void FitShowerAxis(EdbSegP *axis, TObjArray *segments, TCanvas *c1=NULL){
    printf("// Refit axis\n");
    c1->Clear();
    c1->Divide(2,2);
    TGraph grx,gry;
    for(int i=0; i<segments->GetEntriesFast(); i++){
        EdbSegP *s = (EdbSegP *) segments->At(i);
        grx.SetPoint(grx.GetN(), s->Z(), s->X());
        gry.SetPoint(gry.GetN(), s->Z(), s->Y());
    }

    TF1 fx("fx",Form("[0]*(x-%lf)+%lf", axis->Z(), axis->X()));
    TF1 fy("fy",Form("[0]*(x-%lf)+%lf", axis->Z(), axis->Y()));
    fx.SetLineColor(kBlue);
    fy.SetLineColor(kBlue);
    grx.Fit(&fx);
    gry.Fit(&fy);

    
    c1->cd(1);
    grx.Draw("ap");
    float xmin = TMath::MinElement(grx.GetN(), grx.GetX());
    float xmax = TMath::MaxElement(grx.GetN(), grx.GetX());
    TF1 fxo("fxo",Form("%lf*(x-%lf)+%lf", axis->TX(), axis->Z(), axis->X()), xmin, xmax);
    fxo.Draw("same");
    c1->cd(2);
    gry.Draw("ap");
    TF1 fyo("fyo",Form("%lf*(x-%lf)+%lf", axis->TY(), axis->Z(), axis->Y()), xmin, xmax);
    fyo.Draw("same");



    c1->Print("shower.pdf");
    double Z0 = axis->Z();
    axis->SetTX(fx.GetParameter(0));
    axis->SetTY(fy.GetParameter(0));
    axis->PrintNice();

}

void MaskBeamSegments(EdbPVRec *pvr, char *filename, double threshould=250){
    // Flag high density part 
    TFile f(filename, "read");
    TH2F *hTXTY = (TH2F *) f.Get("hTXTY");
    if(hTXTY){
        for(int ipid=0; ipid<pvr->Npatterns();ipid++){
            EdbPattern *pat = pvr->GetPattern(ipid);
            for(int iseg=0; iseg<pat->N(); iseg++){
                EdbSegP *s = pat->GetSegment(iseg);
                double density = hTXTY->GetBinContent(hTXTY->FindBin(s->TX(),s->TY()));
                s->SetFlag( density>threshould?-99:0); // if segment density is higher than 250 segments/cm^2/mrad^2/film do not use.
            }
        }
    }
}

//TNtuple *nt=0;

int CountCylinder(EdbSegP& axis, EdbPVRec *pvr, TNtuple *ntProfile, TNtuple *ntSegments, TObjArray *segments){

    double d2Max=100*100;
    double dt2Max=0.010*0.010;
    double dminMax=50;
    int sum=0;
    int sumDmin=0;

    float Z0 = axis.Z(); // remember the first Z;

    for(int ipid=0; ipid<pvr->Npatterns();ipid++){
        EdbPattern *pat = pvr->GetPattern(ipid);
        float z = pat->Z();
        if(pat->Z()<Z0) continue; 
        if(pat->Z()>Z0+80000) continue; 
        axis.PropagateTo(z);
        int cnt=0;
        int cntDmin=0;
        for(int iseg=0; iseg<pat->N(); iseg++){
            EdbSegP *s = pat->GetSegment(iseg);
            if(s->Flag()==-99) continue; // masked.
            
            // if the segment is upstream of the first Z, continue;
            
            // cylinder and angular cuts
            float d2 = TMath::Power(s->X()-axis.X(),2) + TMath::Power(s->Y()-axis.Y(),2);
            if(d2>d2Max) continue;
            float dt2 = TMath::Power(s->TX()-axis.TX(),2) + TMath::Power(s->TY()-axis.TY(),2);
            if(dt2>dt2Max) continue;

            // dmin cut
            double dminZ;
            float dmin = EdbEDAUtil::CalcDmin(s, &axis, &dminZ);
            if(dminZ<0) dmin=sqrt(d2); // if minz is downstream of the segment, use d.
            
            if(dmin>dminMax) continue;

            cnt++;
            segments->Add(s);
            // if satisfy the cuts, store valuable (if necessary)
            if(ntSegments) ntSegments->Fill(sqrt(d2), sqrt(dt2), dmin);
        }
        sum+=cnt;
        ntProfile->Fill(pat->PID(), cnt);
        printf("ipat %3d, %4d seg, sum %4d seg\n", ipid, cnt, sum);
        
    }
    axis.PropagateTo(Z0);
}

void PlotProfile(TNtuple *ntProfile, TNtuple *ntSegments, TCanvas *c1){

    c1->Clear();
    c1->Divide(2,2);
    int ic1=1;
    c1->cd(ic1++);
    TH1I *hd = new TH1I ("hd","",20,0,100);
    ntSegments->Draw("d >>hd");
    hd->SetXTitle("Position distance from shower axis (#mum)");
    c1->cd(ic1++);
    TH1I *hdt = new TH1I ("hdt","",30,0,0.01);
    ntSegments->Draw("dt >>hdt");
    hdt->SetXTitle("Angular distance from shower axis (rad)");
    c1->cd(ic1++);
    TH1I *hdmin = new TH1I ("hdmin","",20,0,50);
    ntSegments->Draw("dmin >>hdmin");
    hdmin->SetXTitle("Minimum distance from shower axis (#mum)");
    c1->cd(ic1++);
    
    int iplmin = ntProfile->GetMinimum("plate");
    int iplmax = ntProfile->GetMaximum("plate");
    TH1D *hProfile = new TH1D("hProfile", "", (iplmax-iplmin)+1, iplmin-0.5, iplmax+0.5);
    hProfile->SetXTitle("PID (Pattern ID)");
    hProfile->SetYTitle("Number of segments");
    ntProfile->Draw("plate >>hProfile","count*1.0", "goff");
    hProfile->SetMarkerStyle(20);
    hProfile->SetMarkerSize(0.5);
    hProfile->SetLineColor(kBlack);

    // Blur profile to reduce statistical fluctuation
    TH1D *hBlurred = (TH1D *) hProfile->Clone("hBlurred");
    for(int ibin=1; ibin<=hProfile->GetNbinsX(); ibin++){
        double sum = 0;
        for(int i=1; i<=hProfile->GetNbinsX(); i++){
            double val = hProfile->GetBinContent(i);
            sum += val*TMath::Gaus(i,ibin,3.5,1);
        }
        hBlurred->SetBinContent(ibin, sum);
        hBlurred->SetBinError(ibin, 0);
        hProfile->SetBinError(ibin, 0);
    }
    hBlurred->SetLineColor(kBlue);

    hProfile->Draw("p");
    hBlurred->Draw("l, same");
    int maxbin = hBlurred->GetMaximumBin();
    double maxval = hBlurred->GetBinContent(maxbin);
    double maxpid = hBlurred->GetBinCenter(maxbin);
    
    TLine l;
    l.SetLineColor(kBlue);
    l.SetLineWidth(3);
    l.DrawLine(maxpid,0,maxpid, maxval);
    TText tx;
    double showermaximum = maxpid-iplmin+1;
    tx.DrawText(maxpid,maxval*0.3, Form("shower max at %d", (int) showermaximum));
    
    printf("maxbin = %d, maxval = %.1lf, maxpid = %.1lf deltaPID = %.1lf\n", maxbin, maxval, maxpid, showermaximum);
    shower.ShowerMaximum = showermaximum;
    shower.NsegShowerMax = maxval;
    //shower.Erec = ErecShowerMaximum(showermaximum);
    
    double sum = 0, sumShowerMax=0;
    for(int ibin=1; ibin<=hProfile->GetNbinsX(); ibin++){
        double val = hProfile->GetBinContent(ibin);
        sum+=val;
        if(maxbin-3<=ibin&&ibin<=maxbin+3) sumShowerMax+=val;
    }
    shower.NsegSum=sum;
    shower.NsegSumShowerMax=sumShowerMax;

    shower.Erec = ErecNsegSumShowerMax(sumShowerMax);
        
  

}
