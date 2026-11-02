// Analyse final-state particles from the event tree

#include <cmath>
#include <iomanip>
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

// create a histogram that belongs to us, not to the current ROOT file;
// the title has the form "title;x-axis label;y-axis label"
std::unique_ptr<TH1D> makeHistogram(const std::string& name, const std::string& title,
                                    int nBins, double min, double max)
{
    auto histogram = std::make_unique<TH1D>(name.c_str(), title.c_str(), nBins, min, max);
    histogram->SetDirectory(nullptr);
    // keep the sum of squared weights, so bin errors stay correct after scaling
    histogram->Sumw2();
    return histogram;
}

struct CategoryHistograms
{
    std::unique_ptr<TH1D> n;
    std::unique_ptr<TH1D> pt;
    std::unique_ptr<TH1D> eta;
    std::unique_ptr<TH1D> phi;
};

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

// values outside the histogram range are kept in the underflow and overflow bins; report them
void printOutOfRange(const TH1D& histogram)
{
    // every fill has weight 1, so these sums of weights are whole numbers of entries
    const auto entries = static_cast<long long>(histogram.GetEntries());
    const auto underflow = static_cast<long long>(histogram.GetBinContent(0));
    const auto overflow = static_cast<long long>(histogram.GetBinContent(histogram.GetNbinsX() + 1));
    std::cout << std::left << std::setw(10) << histogram.GetName() << std::right
              << " entries " << std::setw(8) << entries
              << "  underflow " << std::setw(7) << underflow
              << " (" << std::setw(6) << 100.0 * underflow / entries << "%)"
              << "  overflow " << std::setw(7) << overflow
              << " (" << std::setw(6) << 100.0 * overflow / entries << "%)\n";
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

    // 4 create four histograms for every category
    std::vector<CategoryHistograms> histograms;
    for (const Category& category : categories)
    {
        histograms.push_back(makeCategoryHistograms(category));
    }

    // 5 read the events one by one and fill the histograms
    long long nEvents = 0;
    long long nParticles = 0;
    {
        // the reader uses the tree in its destructor, so this block ends before the file is closed
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

            // one entry per event and category, also when the multiplicity is zero
            for (std::size_t c = 0; c < categories.size(); ++c)
            {
                histograms[c].n->Fill(multiplicity[c]);
            }
        }

        // stopped for a reason other than the end
        if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd)
        {
            std::cerr << "Error while reading " << argv[1] << " at event " << nEvents << "\n";
            return 1;
        }
    }
    inFile->Close();

    // 6 write the histograms to the output file and close it
    for (const CategoryHistograms& h : histograms)
    {
        outFile->WriteTObject(h.n.get());
        outFile->WriteTObject(h.pt.get());
        outFile->WriteTObject(h.eta.get());
        outFile->WriteTObject(h.phi.get());
    }
    outFile->Close();

    // 7 print the summary
    std::cout << "Events: " << nEvents << ", final-state particles: " << nParticles << "\n";
    std::cout << "Values outside the histogram ranges:\n" << std::fixed << std::setprecision(3);
    for (const CategoryHistograms& h : histograms)
    {
        printOutOfRange(*h.n);
        printOutOfRange(*h.pt);
        printOutOfRange(*h.eta);
        printOutOfRange(*h.phi);
    }
    return 0;
}
