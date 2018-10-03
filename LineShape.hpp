//===============================================
//
//             Master NPAC
// Functions describing Z and J/psi line shape
//
//===============================================


#ifndef __ZLINESHAPE__
#define __ZLINESHAPE__

// Relativistic Brei-Wigner
double BWrel(double *xx, double *par);

// Gaussian
double Gaussian(double *x, double *par);

// Breit Wigner convoluted with a Gaussian
double GBWrel(double *xx, double *par);

double GPapadGBW(double *xx, double *par);

double FitEtaFunc(double* xx, double *par);

#endif
