//===============================================
//
//             Master NPAC
// Implementation of the standard method
// Inherit from FitterBase class
//
//===============================================

#include "TMinuit.h"

#include "FitterStandard.hpp"
#include "TMath.h"

using namespace std;

static FitterStandard* static_pointer_standard = NULL;

//______________________________________________________________________________
void fcn_wrapper_standard(int& npar, double* gin, double& f, double* par, int iflag) {
    if (static_pointer_standard == NULL) {
        cout << "static_pointer_standard is NULL!!!!!" << endl;
        //    exit(-1);
    }
    static_pointer_standard->fcnForMinuit(npar, gin, f, par, iflag);
}


//______________________________________________________________________________
FitterStandard::FitterStandard(MappingTool* map, string name) :
    FitterBase(map, name)
{
}


//______________________________________________________________________________
void FitterStandard::execute()
{
    //==========================================
    //Init
    //==========================================
    static_pointer_standard = this;
    int npar = m_map->getNbOfBins();


    //==========================================
    //TMinuit settings
    //==========================================
    TMinuit minuit(npar);
    minuit.SetPrintLevel(-1);
    minuit.SetFCN(fcn_wrapper_standard);

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
        sprintf(tmp, "%s%d", "alpha", iPar);
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
              cout<<"alpha"<<iPar<<"  "<<m_alphas[iPar] <<"  "<<m_alphaErs[iPar]<<endl;
    }
}


//______________________________________________________________________________
void FitterStandard::fcnForMinuit(int& npar, double* deriv, double& f, double par[], int flag)
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
      //int index0 = m_map->getIndex(m_infoVector->at(i).eta0);
      //int index1 = m_map->getIndex(m_infoVector->at(i).eta1);
      likelihood-= 2*TMath::Log(TMath::Voigt(m_infoVector->at(i).mass/(1+(par[m_infoVector->at(i).region0]+par[m_infoVector->at(i).region1])/2)
       - JPsi_mass, m_infoVector->at(i).sigma, JPsi_width));
    }

    f = likelihood;
}
