//===============================================
//
//             Master NPAC
// Implementation of the BetaAlpha method
// Inherit from FitterBase class
//
//===============================================

#include "TMinuit.h"

#include "FitterBetaAlpha.hpp"
#include "TMath.h"

using namespace std;

static FitterBetaAlpha* static_pointer_BetaAlpha = NULL;

//______________________________________________________________________________
void fcn_wrapper_BetaAlpha(int& npar, double* gin, double& f, double* par, int iflag) {
    if (static_pointer_BetaAlpha == NULL) {
        cout << "static_pointer_BetaAlpha is NULL!!!!!" << endl;
        //    exit(-1);
    }
    static_pointer_BetaAlpha->fcnForMinuit(npar, gin, f, par, iflag);
}


//______________________________________________________________________________
FitterBetaAlpha::FitterBetaAlpha(MappingTool* map, string name) :
    FitterBase(map, name)
{
}


//______________________________________________________________________________
void FitterBetaAlpha::execute()
{
    //==========================================
    //Init
    //==========================================
    static_pointer_BetaAlpha = this;
    int npar = 0;
    int dim = m_map->getNbOfBins();
    //npar=dim*dim
    for(int ivar=dim;ivar>0;ivar--)
    {
      npar+=ivar;
    }
    npar = 576;
    cout << "npar= " << npar << endl;

    //==========================================
    //TMinuit settings
    //==========================================
    TMinuit minuit(npar);
    minuit.SetPrintLevel(-1);
    minuit.SetFCN(fcn_wrapper_BetaAlpha);

    double par[npar];               // the start values
    double stepSize[npar];          // step sizes
    double minVal[npar];            // minimum bound on parameter
    double maxVal[npar];            // maximum bound on parameter
    //string parName[npar];
    string *parName = new string[npar];

    for (int iPar = 0; iPar < npar; iPar++) {
        char tmp[100];
        par[iPar]      = 0;
        stepSize[iPar] = 0.001;
        minVal[iPar]   = -0.3;
        maxVal[iPar]   = +0.3;
        sprintf(tmp, "%s%d%d", "beta_", iPar/24, iPar%24);
        parName[iPar] = tmp;
    }

    for (int iPar=0; iPar<npar; iPar++) {
        minuit.DefineParameter(iPar, parName[iPar].c_str(), par[iPar], stepSize[iPar], minVal[iPar], maxVal[iPar]);
    }

    //==========================================
    //Do the minimization
    //==========================================

    minuit.Migrad();

    //==========================================
    //Get results
    //==========================================

    for (int iPar = 0; iPar < npar; iPar++) {
        minuit.GetParameter(iPar, m_alphas[iPar], m_alphaErs[iPar]);
              cout<<"beta"<<iPar<<"  "<<m_alphas[iPar] <<"  "<<m_alphaErs[iPar]<<endl;
    }
cout<< "&&&&&&&  I'm before constuctAlphas();   &&&&&&&"<<endl;
    constuctAlphas();
cout<< "&&&&&&&  END execute  &&&&&&&"<<endl;
}


//______________________________________________________________________________
void FitterBetaAlpha::fcnForMinuit(int& npar, double* deriv, double& f, double par[], int flag)
{
    //===================================
    //Write here the function to minimize
    //===================================

    //Fonction to minimize
    double likelihood = 0;

    //Definition of the Z-Parameters
    const double Z_mass = 91.1876;
    const double Z_width = 2.4952;

    const double JPsi_mass = 3.096916;
    const double JPsi_width = 0;

    for (unsigned int i = 0; i < m_infoVector->size(); i++)
    {
      int index0 = m_infoVector->at(i).region0;
      int index1 = m_infoVector->at(i).region1;
      //We map the (i,j) doublet on a single parameter
      int index01 = 24 * index0 + index1;
      likelihood-= 2*TMath::Log(TMath::Voigt(m_infoVector->at(i).mass/(1+(par[index01])/2)
       - JPsi_mass, m_infoVector->at(i).sigma, JPsi_width));
    }

    f = likelihood;
}
//______________________________________________________________________________
void FitterBetaAlpha::constuctAlphas()
/*
Reconstruction the alphas values starting from the betas
*/
{


  double betaM[24];
  Float_t U[24*24];

  //get the results
  vector<double> betas_standard   = getAlphas();
  vector<double> betaErs_standard   = getAlphaErs();
  int nBin = m_map->getNbOfBins();
  //We initialize the table to 0
  for(int i=0; i<nBin; i++)
  {
    for(int j=0; j<nBin; j++)
    {
      betaM[i] = 0;
      U[24*i+j] = 0;
    }
  }

  //We fill the table
  for(int i=0; i<nBin; i++)
  {
    for(int j=0; j<nBin; j++)
    {
      double b = betas_standard.at(24*i+j);
      double bEr = betaErs_standard.at(24*i+j);
      if (b != 0)
      {
        betaM[i] += b/(bEr*bEr);
        if(i !=j) {  U[24*i+j]  = 1/(bEr*bEr);}
        else      {  for(int k=0; k < nBin; k ++) {U[24*i+j] += 1/(betaErs_standard.at(24*i+k)*betaErs_standard.at(24*i+k));}}
      }
    }
    //cout << i << endl;
  }
  TMatrix Umat(nBin, nBin, U);
  Umat = Umat.Invert();
  double elem = 0;
  for(int i = 0; i < 24; i++)
  {
    elem = 0;
    for(int j = 0; j < 24; j++)
    {
      elem+= Umat[i][j] * betaM[j];
    }
    alphas.push_back(elem);
    alphaErr.push_back(sqrt(Umat[i][i]));
  }

}
