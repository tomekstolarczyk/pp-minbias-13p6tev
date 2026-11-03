// Generates pp collisions with Pythia 8 and saves the final-state particles to a ROOT tree.
// Usage: generate <settings.cmnd> <output.root>

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"

int main(int argc, char* argv[])
{
    // check that a config file and an output file were given
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <settings.cmnd> <output.root>\n";
        return 1;
    }

    // create the generator object and read the settings file
    Pythia8::Pythia pythia;
    if (!pythia.readFile(argv[1]))
    {
        std::cerr << "Could not read settings from " << argv[1] << "\n";
        return 1;
    }
    if (!pythia.init())
    {
        std::cerr << "Pythia initialization failed\n";
        return 1;
    }

    // create the output file given on the command line and an empty tree "events" inside it
    std::unique_ptr<TFile> outFile(TFile::Open(argv[2], "RECREATE"));
    if (!outFile)
    {
        std::cerr << "Could not create " << argv[2] << "\n";
        return 1;
    }
    auto tree = std::make_unique<TTree>("events", "Final-state particles, one entry per event");

    // connect one list per quantity to a branch; each list holds one value per particle
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

    // generate until the requested number of good events is stored
    const int nEvents = pythia.mode("Main:numberOfEvents");
    const int nAbortedMax = pythia.mode("Main:timesAllowErrors");
    int nAccepted = 0;
    int nAborted = 0;
    long long nFinal = 0;
    while (nAccepted < nEvents)
    {
        // a failed event is inconsistent and must not be used; many failures mean a broken setup
        if (!pythia.next())
        {
            ++nAborted;
            if (nAborted > nAbortedMax)
            {
                std::cerr << "More than " << nAbortedMax << " failed events, stopping the run\n";
                break;
            }
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
            const Pythia8::Particle& particle = pythia.event[i];
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
        tree->Fill();
        ++nAccepted;
    }

    // kOverwrite replaces the copies of the tree header saved automatically while filling
    tree->Write("", TObject::kOverwrite);
    tree.reset();
    outFile->Close();

    pythia.stat();
    std::cout << "Requested events: " << nEvents << ", accepted: " << nAccepted
              << ", aborted: " << nAborted << "\n";
    if (nAccepted > 0)
    {
        std::cout << "Final-state particles: " << nFinal << ", mean per event: " << std::fixed
                  << std::setprecision(2) << static_cast<double>(nFinal) / nAccepted << "\n";
    }
    // a non-zero exit code stops run_all.sh when the run is incomplete
    return nAccepted == nEvents ? 0 : 1;
}
