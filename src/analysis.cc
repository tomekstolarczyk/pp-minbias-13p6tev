// Analyse final-state particles from the event tree

#include <iostream>
#include <memory>
#include <vector>

#include "TFile.h"
#include "TTree.h"

int main(int argc, char* argv[])
{
    // 1 check that an input file was given
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <events.root>\n";
        return 1;
    }

    // 2 open the input file and find the tree "events" inside it
    std::unique_ptr<TFile> inFile(TFile::Open(argv[1], "READ"));
    if (!inFile)
    {
        std::cerr << "Could not open " << argv[1] << "\n";
        return 1;
    }
    // the file owns the tree, so a plain pointer is enough
    TTree* tree = inFile->Get<TTree>("events");
    if (!tree)
    {
        std::cerr << "No tree \"events\" in " << argv[1] << "\n";
        return 1;
    }

    // 3 connect a branch to a vector we own; ROOT fills it on every GetEntry
    std::vector<int> pdg;
    std::vector<int>* pdgAddress = &pdg;
    tree->SetBranchAddress("pdg", &pdgAddress);

    // 4 read the events one by one and count the particles
    const long long nEvents = tree->GetEntries();
    long long nParticles = 0;
    for (long long iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        tree->GetEntry(iEvent);
        nParticles += pdg.size();
    }

    // 5 the tree keeps the address of pdgAddress; forget it before the vector is gone, then close the file
    tree->ResetBranchAddresses();
    inFile->Close();

    // 6 print the summary
    std::cout << "Events: " << nEvents << ", final-state particles: " << nParticles << "\n";
    return 0;
}
