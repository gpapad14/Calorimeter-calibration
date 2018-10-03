#include "LineShape.hpp"


#include <cmath>
#include "TMath.h"
using namespace std;

//===========================================
// Relativistic Brei-Wigner
//===========================================
double BWrel(double *xx, double *par)
/*
  par[0] is the width
  par[1] is the mass
  par[2] is the amplitude
*/
{
   Float_t x =xx[0];
   Double_t f = par[2]*(par[0]*par[0]*x*x/(par[1]*par[1]))/((x*x-par[1]*par[1])*(x*x-par[1]*par[1])+(par[0]*par[0]*x*x*x*x)/(par[1]*par[1]));
   return f;
}

//===========================================
// Gaussian
//===========================================
double Gaussian(double *xx, double *par)
/*
  par[0] is the mean
  par[1] is the sigma
  par[2] is the amplitude
*/
{
     Float_t x =xx[0];
     return par[2]*TMath::Exp(-(x-par[0])*(x-par[0])/(2*par[1]*par[1]));
}

//===========================================
// Breit Wigner convoluted with a Gaussian
//===========================================

double GBWrel(double *xx, double *par)
/*
  par[0] is the sigma of the gaussian function
  par[1] is proportional to the width
  par[2] is the mass
  par[3] is the amplitude
*/
{
   Float_t x =xx[0];
   Double_t f = par[3]*TMath::Voigt(x-par[2], par[0], par[1]);
   return f;
}


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// Breit Wigner convoluted with a Gaussian ... by Giorgos Papadopoulos

double GPapadGBW(double *xx, double *par)
{
   Float_t x =xx[0];
   Double_t f = par[3]*TMath::Pi()*(par[1]/par[2])*(1/sqrt(2*TMath::Pi()) *par[0] ) *
   TMath::Voigt(x*x-par[2]*par[2] , par[0], x*x*2*par[1]/par[2]);
   // [0]: sigma, [1]: gamma, [2]: mass, [3]: amplitude
   return f;
}

double FitEtaFunc(double* xx, double *par)
{
  double x = xx[0];
  return par[0] + par[1]*x + par[2] * x*x;
}
