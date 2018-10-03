//===============================================
//
//             Master NPAC
// Main program for the calibration project
//
//===============================================

// c++ includes
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// ROOT includes
#include "TChain.h"

#include "AnaCalib.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    //command line option
    string inFileName  = "test.list";
    string outFileName = "test.root";
    if(argc>=2) inFileName  = argv[1];
    if(argc>=3) outFileName = argv[2];

    // tree name
    string inTreeName = "tuple";
    // string inTreeName = "tree_ee";

    //=============================================
    //loop over the input files and make a TChain
    //=============================================
    TChain* inChain = new TChain(inTreeName.c_str()); // Need a pointer since AnaCalib ctor is TTree* polymorph.
    cout << "CaloCalib::main()      Load Chain from file: " << inFileName << endl;
    std::ifstream fList(inFileName.c_str());
    if (!fList.is_open()) {
        cerr << "CaloCalib::main()      Can't open file " << inFileName << endl;
        return 1;
    }

    char lineFromFile[10000];
    while (fList.getline(lineFromFile,10000)) {
        TString fileName = lineFromFile;
        if ( inChain->Add(fileName.Data()) ) cout << "CaloCalib::main()      File '" << fileName << "' has been loaded." << endl;
        else cerr << "CaloCalib::main()      Can't load file '" << fileName << "'." << endl;
    }

    cout << "CaloCalib::main()      Total number of entries: " << inChain->GetEntries() << endl << endl;
    fList.close();

    //=============================================
    // Open the output file.
    // This file will store the TH1F, TF1, ... that are created in AnaCalib::Loop().
    //=============================================
    TFile outFile(outFileName.c_str(), "RECREATE");

    //=============================================
    // Do the jobs
    //=============================================
    AnaCalib ana(inChain);
    ana.Loop();

    //=============================================
    //Write in outFile
    //=============================================
    outFile.cd();
    outFile.Write();
    outFile.Close();

    return 0;
}
