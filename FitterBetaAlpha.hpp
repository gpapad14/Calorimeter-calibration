#ifndef FitterBetaAlpha_hpp
#define FitterBetaAlpha_hpp

#include "FitterBase.hpp"
#include "TMatrix.h"
#include<vector>

class FitterBetaAlpha: public FitterBase
{

protected :

    std::vector<double> alphas;
    std::vector<double> alphaErr;

public :

    FitterBetaAlpha(MappingTool* map = NULL, std::string name = "TOTO");

    //Virtual method, to be implemented in daughter class
    virtual void execute();
    virtual void fcnForMinuit(int& npar, double* deriv, double& f, double par[], int flag);
    void constuctAlphas();

    std::vector<double> getAlphaTable() {return alphas;}
    std::vector<double> getAlphaErrTable() {return alphaErr;}
};


#endif
