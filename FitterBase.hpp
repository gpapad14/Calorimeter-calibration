#ifndef FitterBase_hpp
#define FitterBase_hpp

//===============================================
//
//             Master NPAC
// Virtual class common to all fitting methods
// Each fitter inherits from this class
//
// the calibration coefficient are stored in a std::vector
//
//===============================================

// c++ includes
#include <string>

// ROOT inclues
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TGraph.h>

#include "MappingTool.hpp"
#include "InfoForFitter.h"
#include "LineShape.hpp"


class FitterBase
{
public:
    // Constructor.
    FitterBase(MappingTool* map = NULL, std::string name = "TOTO");

    // Destructor.
    virtual ~FitterBase() {};

    // Setter for data.
    void setData(std::vector<InfoForFitter>* infoVector);

    // Get alphas.
    std::vector<double> getAlphas() { return m_alphas; }

    // Get errors on alphas.
    std::vector<double> getAlphaErs() { return m_alphaErs; }

    // Virtual method where the fit is performed, to be implemented in a daughter classes.
    virtual void execute() = 0;

    // Fcn to minimize by Minuit, to be implemented in daughter classes.
    virtual void fcnForMinuit(int& npar, double* deriv, double& f, double par[], int flag) = 0;

    // Write the alphas into a text file.
    void alphaToFile(std::string fileName);


protected :

    // Name associated to the instance of the class.
    std::string m_name;

    // Mapping tool.
    MappingTool* m_map;

    // Fitted values.
    std::vector<double> m_alphas;    // Calibration factors...
    std::vector<double> m_alphaErs;  // ...and their errors.

    // Data.
    std::vector<InfoForFitter>* m_infoVector;

};


#endif
