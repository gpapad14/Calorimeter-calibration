#define AnaCalib_cxx
#include "AnaCalib.hpp"

//root
#include <TLorentzVector.h>
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMatrix.h"

//c++
#include <iostream>
#include <vector>
using namespace std;

//my files
#include "MappingTool.hpp"
#include "FitterStandard.hpp"
#include "FitterBetaAlpha.hpp"
#include "LineShape.hpp"

void AnaCalib::Loop()
{

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast(); // Number of entries -> for instance data_Z_#nentries#

    //data vector
    vector<InfoForFitter> infoVector;

    //Instanciate the mapping tool
    MappingTool map(24, -2.4, 2.4);

    // Other example of instanciation
    vector<double> table = {-2.4, -2., -1.5, -1, -0.8, -0.2, 0.2, 0.8, 1., 1.5, 2., 2.4};
    MappingTool map2(table);

    //Create an histogram for invariant mass
    TH1F* mass_histo = new TH1F("Invariant_Mass", "Mass", 200, 0, 6);

    TH1F* true_mass_histo = new TH1F("Invariant_True_Mass", "True_Mass", 200, 0, 6);
    double alphaFitParam[3]; //Parameters to fit the alpha values at the end of the code

    // Creation of a booleen to decide if we cook the data or not (verification of the code)
    bool cookData = false;
    bool useBeta = true;
    //===============================================
    // Loop over the events
    //===============================================
    //nentries = 3000;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        //Cooking the Energies of electrons
        if(cookData)
        {
          ele0_et*=(1+0.01*pow(-1, map.getIndex(ele0_eta)));
          ele1_et*=(1+0.01*pow(-1, map.getIndex(ele1_eta)));
        }
        //========================================
        //Z events preselection
        //========================================
        //if (ele0_et < 20) continue; //For Z boson
        //if (ele1_et < 20) continue; // For Z boson
        if (fabs(ele0_eta) > 2.4) continue;
        if (fabs(ele1_eta) > 2.4) continue;

        if (ele0_et < 2) continue; // For J/Psi
        if (ele1_et < 2) continue; // For J/Psi
        //========================================
        // Compute the invariant mass: TODO
        //========================================
        double mass = 0;
        mass = sqrt(2*ele0_et*ele1_et*(cosh(ele0_eta-ele1_eta) - cos(ele0_phi-ele1_phi)));
        //if(mass > 80 && mass < 100) mass_histo->Fill(mass); //For Z Boson
        if(mass > 0 && mass < 6) mass_histo->Fill(mass); //For J/Psi

        //========================================
        //Fill InfoForFitter
        //========================================

        InfoForFitter info;

        info.eta0     = ele0_eta;
        info.phi0     = ele0_phi;
        info.pt0      = ele0_et;
        info.region0  = map.getIndex(ele0_eta);
        info.eta1     = ele1_eta;
        info.phi1     = ele1_phi;
        info.pt1      = ele1_et;
        info.region1  = map.getIndex(ele1_eta);
        info.mass     = mass;
        infoVector.push_back(info);

    } //end of the event loop

    //========================================
    //Fit
    //========================================

    //We fit the invariant mass distribution


    //Function to fit (GBWrel) is in LineShape.cpp
    /* Z boson
    TF1* f1 = new TF1("myGBW", GBWrel, 80, 100, 4);
    f1->SetParameters(1, 4, 91, 700); // for mc 10000
    f1->SetParNames("sigma", "width", "mass", "amplitude");
    mass_histo->Fit("myGBW");
    */
    //J/Psi
    TF1* f1 = new TF1("gauss", Gaussian, 0, 6, 3);
    f1->SetParameters(3, 0.1, 5000); // for mc 10000
    f1->SetParNames( "mass", "sigma", "amplitude");
    mass_histo->Fit("gauss");

    /*
    TF1* f1 = new TF1("fitfunc", GPapadGBW, 80, 100, 4);
    f1->SetParameters(1, 2, 91, 1);
    f1->SetParNames("sigma", "width", "mass", "amplitude");
    mass_histo->Fit("fitfunc");
    */

    //TCanvas* cf2 = new TCanvas("cf2","GBW function");
    //TF1* f2 = new TF1("truth_myGBW", GBWrel, 80, 100, 4);
    //f2->SetParameters(1, 4, 91, 700);
    //f2->SetParNames("sigma", "width", "mass", "amplitude");
    //f2->Draw("");
    //cf2->SaveAs("GBWdistribution.root");
    //truth_mass_histo->Fit("truth_myGBW");

    //We extract the value of the sigma parameter
    double sigma = f1->GetParameter("sigma");
    //cout << sigma << endl;

    for(unsigned int i=0; i<infoVector.size(); i++)
      {
        infoVector.at(i).sigma = sigma;
      }

    vector<double> alphas_standard;
    vector<double> alphaErs_standard;
    vector<double> betas_standard;
    vector<double> betaErs_standard;

    if(!useBeta)
    /*
    If we the direct alpha parametrization Method
    */
    {
      //initialize the fitter tool
      FitterStandard fit_standard(&map, "std");
      fit_standard.setData(&infoVector);
      // do the fit
      fit_standard.execute();

      //get the results
      alphas_standard   = fit_standard.getAlphas();
      alphaErs_standard = fit_standard.getAlphaErs();
    }
    else
    /*
    If we the beta parametrization Method
    */
    {
      //initialize the fitter tool
      FitterBetaAlpha fit_BetaAlpha(&map, "std");
      fit_BetaAlpha.setData(&infoVector);
      // do the fit
cout<< "*********    I'm before fit_BetaAlpha.execute()    *********"<<endl;
      fit_BetaAlpha.execute();
cout<< "*********    I'm after fit_BetaAlpha.execute()    *********"<<endl;

      alphas_standard = fit_BetaAlpha.getAlphaTable();
      alphaErs_standard = fit_BetaAlpha.getAlphaErrTable();
      /*
      for(unsigned int i = 0; i<alphas_standard.size(); i++)
      {
        cout << alphas_standard.at(i) << endl;
      }
      */
    }

    TCanvas* c1 = new TCanvas();

    //Defining all the table that could be used in the rest of the code
    vector<double>* index = new vector<double>;
    vector<double>* indexErr = new vector<double>;
    vector<double>* diff_value = new vector<double>;
    vector<double>* pull_value = new vector<double>;
    vector<double>* alpha_value = new vector<double>;
    vector<double>* alphaErr_value = new vector<double>;


    for(int i = 0; i<(int)(alphas_standard.size()); i++)
    {
      double eta = 0;
      map.getEta(i, eta);
      index->push_back(eta);
      indexErr->push_back(0.1);
      if(cookData)
      /*
      If we are using MC simulation data and adding an error in the data
      */
      {
        diff_value->push_back(alphas_standard.at(i)-0.01*pow(-1, i));
        pull_value->push_back((alphas_standard.at(i)-0.01*pow(-1, i))/alphaErs_standard.at(i));
      }
      else
      /*
      If we are using real data
      */
      {
        alpha_value->push_back(alphas_standard.at(i));
        alphaErr_value->push_back(alphaErs_standard.at(i));
      }
    }

    if(cookData)
    /*
    If we are using MC simulation data and adding an error in the data
    */
    {
      c1 -> Divide(2,1);
      c1->cd(1);
      TGraph* diff = new TGraph((Int_t)(index->size()), index->data(), diff_value->data());
      diff->SetTitle("a_(i) - lambda");
      diff->Draw("AP*");
      c1->cd(2);
      TGraph* pull = new TGraph((Int_t)(index->size()), index->data(), pull_value->data());
      pull->Draw("AP*");
      c1->SaveAs("Results_Bias.root");
    }
    else
    /*
    If we are using real data
    */
    {

      TF1* fiteta = new TF1("polfeta2", FitEtaFunc, -2.4, 2.4, 3);
      fiteta->SetParameters(-0.01, -0.0001, -0.01); // for mc 10000
      fiteta->SetParNames("p0", "p1", "p2");

      TGraphErrors* alpha = new TGraphErrors((Int_t)(index->size()), index->data(), alpha_value->data(), indexErr->data(), alphaErr_value->data());
      alpha->SetTitle("alpha_i");
      alpha->Draw("AP*");
      alpha->Fit("polfeta2");
      c1->SaveAs("Results_JPsi_beta.root");

      alphaFitParam[0] = fiteta->GetParameter(0);
      alphaFitParam[1] = fiteta->GetParameter(1);
      alphaFitParam[2] = fiteta->GetParameter(2);


    }

// Filling of the true_mass histogram
      for(unsigned int i = 0; i < infoVector.size(); i++)
      {
        double true_mass=0;

        double alpha0 = 0, alpha1 = 0;
        alpha0 = FitEtaFunc(&(infoVector.at(i).eta0), alphaFitParam);
        alpha1 = FitEtaFunc(&(infoVector.at(i).eta1), alphaFitParam);

        true_mass = infoVector.at(i).mass/(1+(alpha0 + alpha1)/2);
        true_mass_histo->Fill(true_mass);
      }

/* Z-boson
// Fitting of the true_mass (calibrated) by GBW function
      TF1* true_f1 = new TF1("myGBW_true", GBWrel, 80, 100, 4);
      true_f1->SetParameters(1, 4, 91, 700); // for mc 10000
      true_f1->SetParNames("sigma", "width", "mass", "amplitude");
      true_mass_histo->Fit("myGBW_true");
      */
/* JPsi */
    // Fitting of the true_mass (calibrated) by GBW function
      TF1* true_f1 = new TF1("gauss_true", Gaussian, 0, 6, 3);
      true_f1->SetParameters(3, 0.1, 5000); // for mc 10000
      true_f1->SetParNames( "mass", "sigma", "amplitude");
      true_mass_histo->Fit("gauss_true");

} //end of AnaCalib::Loop()
