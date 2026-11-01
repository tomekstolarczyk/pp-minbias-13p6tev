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
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        if (!pythia.next())
        {
            ++nAborted;
        }
    }

    // 5 print the summary
    pythia.stat();
    std::cout << "Requested events: " << nEvents << ", aborted: " << nAborted << "\n";
    return 0;
}
