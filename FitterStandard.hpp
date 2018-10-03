#ifndef FitterStandard_hpp
#define FitterStandard_hpp

#include "FitterBase.hpp"

class FitterStandard: public FitterBase
{

protected :

public :

    FitterStandard(MappingTool* map = NULL, std::string name = "TOTO");

    //Virtual method, to be implemented in daughter class
    virtual void execute();
    virtual void fcnForMinuit(int& npar, double* deriv, double& f, double par[], int flag);

};


#endif
