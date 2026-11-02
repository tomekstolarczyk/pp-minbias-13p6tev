// Generate pp collisions with Pythia 8

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
using namespace Pythia8;

int main(int argc, char* argv[])
{
    // 1 check that a settings file and an output file were given
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <settings.cmnd> <output.root>\n";
        return 1;
    }

    // 2 create the generator object and read the settings file
    Pythia pythia;
    if (!pythia.readFile(argv[1]))
    {
        std::cerr << "Could not read settings from " << argv[1] << "\n";
        return 1;
    }

    // 3 prepare the generator
    if (!pythia.init())
    {
        std::cerr << "Pythia initialization failed\n";
        return 1;
    }

    // 4 create the output file given on the command line and an empty tree "events" inside it
    std::unique_ptr<TFile> outFile(TFile::Open(argv[2], "RECREATE"));
    if (!outFile)
    {
        std::cerr << "Could not create " << argv[2] << "\n";
        return 1;
    }
    auto tree = std::make_unique<TTree>("events", "Final-state particles, one entry per event");

    // 5 connect one list per quantity to a branch; each list holds one value per particle
    std::vector<int> pdg;
    std::vector<int> charge;
    std::vector<double> pt;
    std::vector<double> eta;
    std::vector<double> phi;
    tree->Branch("pdg", &pdg);
    tree->Branch("charge", &charge);
    tree->Branch("pt", &pt);
    tree->Branch("eta", &eta);
    tree->Branch("phi", &phi);

    // 6 generate the events
    const int nEvents = pythia.mode("Main:numberOfEvents");
    int nAborted = 0;
    long long nFinal = 0;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        // an aborted event is incomplete, so it must not be used
        if (!pythia.next())
        {
            ++nAborted;
            continue;
        }

        pdg.clear();
        charge.clear();
        pt.clear();
        eta.clear();
        phi.clear();

        // the event record also holds intermediate particles, keep only the final ones
        for (int i = 0; i < pythia.event.size(); ++i)
        {
            const Particle& particle = pythia.event[i];
            if (particle.isFinal())
            {
                ++nFinal;
                pdg.push_back(particle.id());
                // chargeType() is three times the charge; final-state particles have integer charge
                charge.push_back(particle.chargeType() / 3);
                pt.push_back(particle.pT());
                eta.push_back(particle.eta());
                phi.push_back(particle.phi());
            }
        }

        // store the lists as one entry of the tree, also when they are empty
        tree->Fill();
    }

    // 7 save the tree to file, delete the tree and then close the file
    // kOverwrite replaces the header copies that the tree saved automatically while filling
    tree->Write("", TObject::kOverwrite);
    tree.reset();
    outFile->Close();

    // 8 print the summary
    pythia.stat();
    const int nAccepted = nEvents - nAborted;
    std::cout << "Requested events: " << nEvents << ", aborted: " << nAborted << "\n";
    std::cout << "Final-state particles: " << nFinal
              << ", per event: " << std::fixed << std::setprecision(2)
              << static_cast<double>(nFinal) / nAccepted << "\n";
    return 0;
}
