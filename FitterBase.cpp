#include "FitterBase.hpp"

#include <fstream>
#include <iomanip>

using namespace std;

//______________________________________________________________________________
FitterBase::FitterBase(MappingTool* map, std::string name)
{
    m_map = map;
    m_name = name;
    m_infoVector = NULL;

    if (map != NULL) {
        //m_alphas.resize( m_map->getNbOfBins() );
        //m_alphaErs.resize( m_map->getNbOfBins() );
        m_alphas.resize( 576 );
        m_alphaErs.resize( 576 );
        for (size_t iAlpha = 0; iAlpha < m_alphas.size(); iAlpha++) {
            m_alphas[iAlpha]     = -999;
            m_alphaErs[iAlpha]   = -999;
        }
    }
}


//______________________________________________________________________________
void FitterBase::setData(std::vector<InfoForFitter>* infoVector)
{
    m_infoVector = infoVector;
}


//______________________________________________________________________________
void FitterBase::alphaToFile(string filename)
{
    std::ofstream outputFile(filename.c_str());

    for (size_t iAlpha = 0; iAlpha < m_alphas.size(); iAlpha++) {
        outputFile << iAlpha << " "
                   << setw(10) << m_alphas[iAlpha]   << "  "
                   << setw(10) << m_alphaErs[iAlpha] << "  "
                   << endl;
    }

    outputFile.close();
}
