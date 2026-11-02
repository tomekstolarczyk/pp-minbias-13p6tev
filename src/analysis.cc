// Analyse final-state particles from the event tree

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TTree.h"

namespace
{

// histogram binning from the task specification
constexpr int nBinsN = 400;
constexpr double minN = 0.0;
constexpr double maxN = 400.0;

constexpr int nBinsPt = 100;
constexpr double minPt = 0.0; // GeV
constexpr double maxPt = 10.0; // GeV

constexpr int nBinsEta = 100;
constexpr double minEta = -5.0;
constexpr double maxEta = 5.0;

constexpr int nBinsPhi = 64;
constexpr double minPhi = 0.0;
constexpr double maxPhi = TMath::Pi();

// create a histogram that belongs to us, not to the current ROOT file
std::unique_ptr<TH1D> makeHistogram(const std::string& name, const std::string& title, int nBins, double min, double max)
{
    auto histogram = std::make_unique<TH1D>(name.c_str(), title.c_str(), nBins, min, max);
    histogram->SetDirectory(nullptr);
    // keep the sum of squared weights, so bin errors stay correct after scaling
    histogram->Sumw2();
    return histogram;
}

}

int main(int argc, char* argv[])
{
    // 1 check that an input file and an output file were given
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <events.root> <analysis.root>\n";
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

    // 3 create the output file now, so a wrong path fails before the long loop
    std::unique_ptr<TFile> outFile(TFile::Open(argv[2], "RECREATE"));
    if (!outFile)
    {
        std::cerr << "Could not create " << argv[2] << "\n";
        return 1;
    }

    // 4 connect the branches to vectors we own; ROOT fills them on every GetEntry
    std::vector<int> pdg;
    std::vector<double> pt;
    std::vector<double> eta;
    std::vector<double> phi;
    std::vector<int>* pdgAddress = &pdg;
    std::vector<double>* ptAddress = &pt;
    std::vector<double>* etaAddress = &eta;
    std::vector<double>* phiAddress = &phi;
    tree->SetBranchAddress("pdg", &pdgAddress);
    tree->SetBranchAddress("pt", &ptAddress);
    tree->SetBranchAddress("eta", &etaAddress);
    tree->SetBranchAddress("phi", &phiAddress);

    // 5 create the histograms of all final-state particles
    auto hN = makeHistogram("hN_all", "Multiplicity of all final-state particles", nBinsN, minN, maxN);
    auto hPt = makeHistogram("hPt_all", "p_{T} of all final-state particles", nBinsPt, minPt, maxPt);
    auto hEta = makeHistogram("hEta_all", "#eta of all final-state particles", nBinsEta, minEta, maxEta);
    auto hPhi = makeHistogram("hPhi_all", "|#phi| of all final-state particles", nBinsPhi, minPhi, maxPhi);

    // 6 read the events one by one and fill the histograms
    const long long nEvents = tree->GetEntries();
    long long nParticles = 0;
    for (long long iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        tree->GetEntry(iEvent);
        nParticles += pdg.size();

        // one entry per event, also when the event has no particles
        hN->Fill(pdg.size());
        for (std::size_t i = 0; i < pdg.size(); ++i)
        {
            hPt->Fill(pt[i]);
            hEta->Fill(eta[i]);
            // phi is in (-pi, pi] but the specified range is [0, pi]; isotropy makes |phi| equivalent
            hPhi->Fill(std::abs(phi[i]));
        }
    }

    // 7 the tree keeps the addresses of our pointers; forget them before the vectors are gone
    tree->ResetBranchAddresses();
    inFile->Close();

    // 8 write the histograms to the output file and close it
    outFile->WriteTObject(hN.get());
    outFile->WriteTObject(hPt.get());
    outFile->WriteTObject(hEta.get());
    outFile->WriteTObject(hPhi.get());
    outFile->Close();

    // 9 print the summary
    std::cout << "Events: " << nEvents << ", final-state particles: " << nParticles << "\n";
    return 0;
}
