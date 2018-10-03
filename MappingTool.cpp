#include "MappingTool.hpp"
#include <vector>
#include <cmath>

using namespace std;

MappingTool::MappingTool()
{
    setEtaBins(24, -2.4, 2.4);
}

MappingTool::MappingTool(const vector<double>& etaBorders)
{
    setEtaBins(etaBorders);
}


MappingTool::MappingTool(int nEtaBins, double etamin, double etamax)
{
    setEtaBins(nEtaBins, etamin, etamax);
}



void MappingTool::setEtaBins(const vector<double>& etaBorders)
{
    m_EtaBorders.clear();
    m_EtaBorders = etaBorders;
}

void MappingTool::setEtaBins(int nEtaBins, double etamin, double etamax)
{
    m_EtaBorders.clear();
    m_EtaBorders.resize(nEtaBins + 1);

    double step = (etamax - etamin) / nEtaBins;
    for (int i = 0; i < nEtaBins + 1; i++)
        m_EtaBorders[i] = etamin + i*step;
}


void MappingTool::print()
{
    cout << "Eta bins = { ";
    for (int i = 0; i < (int)m_EtaBorders.size(); i++) {
        cout << m_EtaBorders[i] << " ";
    }
    cout << "}" << endl;
}

int MappingTool::getIndex(double eta)
{
    int foundEtaBin = -1;
    for (int iEta = 0; iEta < (int)m_EtaBorders.size()-1; iEta++) {
        if (m_EtaBorders[iEta] <= eta && eta < m_EtaBorders[iEta+1]) {
            foundEtaBin = iEta;
            break;
        }
    }

    if (foundEtaBin < 0 || foundEtaBin >= (int)getNbOfBins()) {
        cout << "MappingTool::getIndex()     Could not map eta = " << eta
             << " to a map index. Returning -1." << endl;
    }

    return foundEtaBin;
}


void MappingTool::getEta(int index, double& etaCent)
{
    etaCent = m_EtaBorders[index] + (m_EtaBorders[index+1] - m_EtaBorders[index])/2;
}



double MappingTool::getEtaBinSize(int index)
{
    return m_EtaBorders[index+1] - m_EtaBorders[index];
}
