// Generate pp collisions with Pythia 8

#include <iostream>
#include <memory>

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

    // 5 generate the events
    const int nEvents = pythia.mode("Main:numberOfEvents");
    int nAborted = 0;
    long nFinal = 0;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        // an aborted event is incomplete, so it must not be used
        if (!pythia.next())
        {
            ++nAborted;
            continue;
        }

        // the event record also holds intermediate particles, keep only the final ones
        for (int i = 0; i < pythia.event.size(); ++i)
        {
            if (pythia.event[i].isFinal())
            {
                ++nFinal;
            }
        }
    }

    // 6 save the tree to file, delete the tree and then close the file
    tree->Write();
    tree.reset();
    outFile->Close();

    // 7 print the summary
    pythia.stat();
    const int nAccepted = nEvents - nAborted;
    std::cout << "Requested events: " << nEvents << ", aborted: " << nAborted << "\n";
    std::cout << "Final-state particles: " << nFinal
              << ", per event: " << static_cast<double>(nFinal) / nAccepted << "\n";
    return 0;
}
