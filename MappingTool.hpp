#ifndef MAPPINGTOOL_HPP
#define MAPPINGTOOL_HPP

//===============================================
//
//                 Master NPAC
//
// For each eta bins, an integer is associated and vice-versa
//        (eta)    <--->  integer
//
//
//  Constructor 1:
//      MappingTool map(  24, -2.4,2.4);
//  Constructor 2:
//      vector<double> table;
//      MappingTool map2(table);
//      this constructor allows bins with different size
//
//===============================================

#include <iostream>
#include <vector>

class MappingTool
{
public:

    /** constructor */
    MappingTool();

    /** constructor */
    MappingTool(const std::vector<double>& etaBorders);

    /** constructor */
    MappingTool(int nEta, double etamax, double etamin);

    /** Set Eta bins */
    void setEtaBins(const std::vector<double>& etaBorders);

    /** Set Eta bins */
    void setEtaBins(int nEta, double etamax, double etamin);



    /** get index from (eta) */
    int getIndex(double eta);

    /** get (eta) from index */
    void getEta(int index, double& etaCent);

    /** get the total number of regions */
    unsigned int getNbOfBins() { return m_EtaBorders.size() - 1; }

    /** get the size of the regions in eta */
    double getEtaBinSize(int binIndex);

    /** clear the vectors */
    void clear() { m_EtaBorders.clear(); }

    /** Print the vector. */
    void print();

    /** Get bins. */
    std::vector<double> getEtaBorders() { return m_EtaBorders; }


protected:

    /** Bin borders. */
    std::vector<double> m_EtaBorders;

};

#endif
