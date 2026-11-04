// Fills multiplicity, pT, eta and phi histograms of final-state particles for eight categories.
// Usage: analysis <events.root> <analysis.root>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

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

// PDG codes of the particles; antiparticles have the negative code
constexpr int pdgPion = 211;
constexpr int pdgKaon = 321;
constexpr int pdgProton = 2212;

// decides from the PDG code and the charge whether a particle belongs to a category
using Selector = bool (*)(int pdg, int charge);

struct Category
{
    std::string name;  // used in histogram names, e.g. hPt_pip
    std::string label; // used in histogram titles
    Selector selects;
};

// "all" includes neutral particles, "ch" only charged ones
const std::vector<Category> categories = {
    {"all", "all final-state particles", [](int, int) { return true; }},
    {"ch", "charged particles", [](int, int charge) { return charge != 0; }},
    {"pip", "#pi^{+}", [](int pdg, int) { return pdg == pdgPion; }},
    {"pim", "#pi^{-}", [](int pdg, int) { return pdg == -pdgPion; }},
    {"Kp", "K^{+}", [](int pdg, int) { return pdg == pdgKaon; }},
    {"Km", "K^{-}", [](int pdg, int) { return pdg == -pdgKaon; }},
    {"p", "p", [](int pdg, int) { return pdg == pdgProton; }},
    {"pbar", "#bar{p}", [](int pdg, int) { return pdg == -pdgProton; }},
};

// creates a histogram owned by the caller, not by the current ROOT file;
// the title has the form "title;x-axis label;y-axis label"
std::unique_ptr<TH1D> makeHistogram(const std::string& name, const std::string& title,
                                    int nBins, double min, double max)
{
    auto histogram = std::make_unique<TH1D>(name.c_str(), title.c_str(), nBins, min, max);
    // detach it from the current ROOT file, so that only the unique_ptr deletes it
    histogram->SetDirectory(nullptr);
    // keep the sum of squared weights, so bin errors stay correct after scaling
    histogram->Sumw2();
    // include under- and overflow in the mean and standard deviation, so they are the true values
    histogram->SetStatOverflows(TH1::kConsider);
    return histogram;
}

struct CategoryHistograms
{
    std::unique_ptr<TH1D> n;
    std::unique_ptr<TH1D> pt;
    std::unique_ptr<TH1D> eta;
    std::unique_ptr<TH1D> phi;
};

// four histograms (N, pT, eta, phi) for each category
CategoryHistograms makeCategoryHistograms(const Category& category)
{
    const std::string& name = category.name;
    const std::string& label = category.label;
    return {
        makeHistogram("hN_" + name, "Multiplicity of " + label + ";Multiplicity N;Number of events",
                      nBinsN, minN, maxN),
        makeHistogram("hPt_" + name, "p_{T} of " + label + ";p_{T} [GeV];Number of particles",
                      nBinsPt, minPt, maxPt),
        makeHistogram("hEta_" + name, "#eta of " + label + ";#eta;Number of particles",
                      nBinsEta, minEta, maxEta),
        makeHistogram("hPhi_" + name, "#phi of " + label + ";#phi [rad];Number of particles",
                      nBinsPhi, minPhi, maxPhi),
    };
}

// one line per histogram: the true mean and standard deviation (over all entries, see makeHistogram)
// and the values outside the range, which are kept in the underflow and overflow bins
void printStatistics(const TH1D& histogram)
{
    // ROOT returns these counts as double
    const double entries = histogram.GetEntries();
    const double underflow = histogram.GetBinContent(0);
    const double overflow = histogram.GetBinContent(histogram.GetNbinsX() + 1);
    std::printf("%-10s %9.0f %9.4f %9.4f %8.0f (%6.3f%%) %8.0f (%6.3f%%)\n", histogram.GetName(),
                entries, histogram.GetMean(), histogram.GetStdDev(), underflow,
                100.0 * underflow / entries, overflow, 100.0 * overflow / entries);
}
}

int main(int argc, char* argv[])
{
    // check that an input file and an output file were given
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <events.root> <analysis.root>\n";
        return 1;
    }

    // open the input file and find the tree "events" inside it
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

    // create the output file before the long loop, so that a wrong path fails early
    std::unique_ptr<TFile> outFile(TFile::Open(argv[2], "RECREATE"));
    if (!outFile)
    {
        std::cerr << "Could not create " << argv[2] << "\n";
        return 1;
    }

    // create four histograms for every category
    std::vector<CategoryHistograms> histograms;
    for (const Category& category : categories)
    {
        histograms.push_back(makeCategoryHistograms(category));
    }

    // read the events one by one and fill the histograms
    long long nEvents = 0;
    long long nParticles = 0;
    {
        // the reader uses the tree in its destructor, so this block ends before the input file is closed
        TTreeReader reader(tree);
        // each array holds one value per particle of the current event
        TTreeReaderArray<int> pdg(reader, "pdg");
        TTreeReaderArray<int> charge(reader, "charge");
        TTreeReaderArray<double> pt(reader, "pt");
        TTreeReaderArray<double> eta(reader, "eta");
        TTreeReaderArray<double> phi(reader, "phi");

        while (reader.Next())
        {
            ++nEvents;
            nParticles += pdg.GetSize();

            // number of particles of each category in this event
            std::vector<int> multiplicity(categories.size(), 0);
            for (std::size_t i = 0; i < pdg.GetSize(); ++i)
            {
                // a particle can be in several categories (a pi+ is in all, ch and pip), so check them all
                for (std::size_t c = 0; c < categories.size(); ++c)
                {
                    if (!categories[c].selects(pdg[i], charge[i]))
                    {
                        continue;
                    }
                    ++multiplicity[c];
                    histograms[c].pt->Fill(pt[i]);
                    histograms[c].eta->Fill(eta[i]);
                    // phi is in (-pi, pi] but the specified range is [0, pi]; isotropy makes |phi| equivalent
                    histograms[c].phi->Fill(std::abs(phi[i]));
                }
            }
            // fill the multiplicities once per event, also when they are zero
            for (std::size_t c = 0; c < categories.size(); ++c)
            {
                histograms[c].n->Fill(multiplicity[c]);
            }
        }

        // Next() also returns false on a read error, not only at the end of the tree
        if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd)
        {
            std::cerr << "Error while reading " << argv[1] << " at event " << nEvents << "\n";
            return 1;
        }
    }
    inFile->Close();

    // write the histograms to the output file and close it
    for (const CategoryHistograms& h : histograms)
    {
        outFile->WriteTObject(h.n.get());
        outFile->WriteTObject(h.pt.get());
        outFile->WriteTObject(h.eta.get());
        outFile->WriteTObject(h.phi.get());
    }
    outFile->Close();

    // the particle count can be compared with the one printed by the generator
    std::cout << "Events: " << nEvents << ", final-state particles: " << nParticles << "\n";
    std::cout << "Statistics over all entries, including those outside the histogram range:\n";
    std::printf("%-10s %9s %9s %9s %19s %19s\n", "histogram", "entries", "mean", "std dev",
                "underflow", "overflow");
    for (const CategoryHistograms& h : histograms)
    {
        printStatistics(*h.n);
        printStatistics(*h.pt);
        printStatistics(*h.eta);
        printStatistics(*h.phi);
    }
    return 0;
}
