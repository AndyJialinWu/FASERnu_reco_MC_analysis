#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
#include <TVector3.h>
 
#include <cassert>
 
using namespace ROOT::Math;
 
 
// define the parametric line equation
void line(double t, const double *p, double &x, double &y, double &z) {
   // a parametric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
};
 
 

 
// function Object to be minimized
struct SumDistance2 {
   // the TGraph is a data member of the object
   TGraph2D *fGraph;
   bool first = true;
 
   SumDistance2(TGraph2D *g) : fGraph(g) {}
 
   // calculate distance line-point
   double distance2(double x, double y, double z, const double *p) {
      // distance line point is D= | (xp-x0) cross  ux |
      // where ux is direction of line and x0 is a point in the line (like t = 0)
      XYZVector xp(x,y,z);
      XYZVector x0(p[0], p[2], 0. );
      XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
      XYZVector u = (x1-x0).Unit();
      double d2 = ((xp-x0).Cross(u)).Mag2();
      return d2;
   }
 
   // implementation of the function to be minimized
   double operator() (const double *par) {
      assert(fGraph != nullptr);
      double * x = fGraph->GetX();
      double * y = fGraph->GetY();
      double * z = fGraph->GetZ();
      int npoints = fGraph->GetN();
      double sum = 0;
      for (int i = 0; i < npoints; ++i) {
         double d = distance2(x[i],y[i],z[i],par);
         sum += d;
      }
      if (first) {
         std::cout << "Total Initial distance square = " << sum << std::endl;
      }
      first = false;
      return sum;
   }
 
};

const double *line3Dfit(std::vector<float> X, std::vector<float> Y, std::vector<float> Z){

   int npoints = Z.size();
   if(npoints < 2){
      std::cout<<"The number of points is not enough for defining a line !"<<std::endl;
   }

   if(npoints > 5) npoints = 5; // fit first 5 segments

   TGraph2D * gr = new TGraph2D();
   // generate graph with the 3d points
   for(int i=0; i<npoints; i++){
      gr->SetPoint(i, X[i], Y[i], Z[i]);
   }

   ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
   
   // make the functor objet
   SumDistance2 sdist(gr);
   ROOT::Math::Functor fcn(sdist,4);
   // set the function and the initial parameter values
   //double pStart[4] = {X[0], 0, Y[0], 0};
   double pStart[4] = {1, 1, 1, 1};
   
   if((Z[1] - Z[0]) != 0){
      pStart[1] = (X[1] - X[0]) / (Z[1] - Z[0]);
      pStart[3] = (Y[1] - Y[0]) / (Z[1] - Z[0]);

      pStart[0] = X[0] - pStart[1] * Z[0];
      pStart[2] = Y[0] - pStart[3] * Z[0];
   }
   

   fitter->SetFCN(fcn,pStart);
   // set step sizes different than default ones (0.0001 times parameter values)
   for (int i = 0; i < 4; ++i) fitter->Config().ParSettings(i).SetStepSize(0.0001);

   bool ok = fitter->FitFCN();
   if (!ok) {
      Error("line3Dfit","Line3D Fit failed");
   }

   const ROOT::Fit::FitResult & result = fitter->Result();
 
   std::cout << "Total final distance square = " << result.MinFcnValue() << std::endl;
   //result.Print(std::cout);

   const double *parFit = result.GetParams(); // p
   return parFit;

};

double GetImpactParameter(TVector3 TrkDirVect, double p0, double p2, float vx_hit, float vy_hit, float vz_hit){
   // p0 and p2 are line 3D fitting parameters

   TVector3 xv(vx_hit, vy_hit, vz_hit);
   TVector3 x0(p0, p2, 0. );
   TVector3 u = TrkDirVect.Unit();

   double IP = ((xv-x0).Cross(u)).Mag();
   return IP;

};


