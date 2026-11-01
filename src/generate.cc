// Generate pp collisions with Pythia 8 

#include <iostream>

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main(int argc, char* argv[])
{
    // 1 check that a settings file was given
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <settings.cmnd>\n";
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

    // 4 generate the events
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

    // 5 print the summary
    pythia.stat();
    const int nAccepted = nEvents - nAborted;
    std::cout << "Requested events: " << nEvents << ", aborted: " << nAborted << "\n";
    std::cout << "Final-state particles: " << nFinal
              << ", per event: " << static_cast<double>(nFinal) / nAccepted << "\n";
    return 0;
}
