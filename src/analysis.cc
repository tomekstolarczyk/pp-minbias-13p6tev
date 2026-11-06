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
#include "TNamed.h"
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

constexpr int pdgPion = 211;
constexpr int pdgKaon = 321;
constexpr int pdgProton = 2212;

using Selector = bool (*)(int pdg, int charge);

struct Category
{
    std::string name;
    std::string label;
    Selector selects;
};

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

void printStatistics(const TH1D& histogram)
{
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
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <events.root> <analysis.root>\n";
        return 1;
    }

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
    // the generator settings are copied to the output, so analysis.root records how it was produced;
    // a copy, because the input file deletes the original when it is closed
    const TNamed* settingsInFile = inFile->Get<TNamed>("generator_settings");
    if (!settingsInFile)
    {
        std::cerr << "No generator_settings in " << argv[1] << "\n";
        return 1;
    }
    TNamed generatorSettings(*settingsInFile);

    std::unique_ptr<TFile> outFile(TFile::Open(argv[2], "RECREATE"));
    if (!outFile)
    {
        std::cerr << "Could not create " << argv[2] << "\n";
        return 1;
    }

    std::vector<CategoryHistograms> histograms;
    for (const Category& category : categories)
    {
        histograms.push_back(makeCategoryHistograms(category));
    }

    long long nEvents = 0;
    long long nParticles = 0;
    // the reader uses the tree in its destructor, so this block ends before the input file is closed
    {
        TTreeReader reader(tree);
        TTreeReaderArray<int> pdg(reader, "pdg");
        TTreeReaderArray<int> charge(reader, "charge");
        TTreeReaderArray<double> pt(reader, "pt");
        TTreeReaderArray<double> eta(reader, "eta");
        TTreeReaderArray<double> phi(reader, "phi");

        while (reader.Next())
        {
            ++nEvents;
            nParticles += pdg.GetSize();

            std::vector<int> multiplicity(categories.size(), 0);
            for (std::size_t i = 0; i < pdg.GetSize(); ++i)
            {
                // a particle can be in several categories, so check them all
                for (std::size_t c = 0; c < categories.size(); ++c)
                {
                    if (!categories[c].selects(pdg[i], charge[i]))
                    {
                        continue;
                    }
                    ++multiplicity[c];
                    histograms[c].pt->Fill(pt[i]);
                    histograms[c].eta->Fill(eta[i]);
                    // phi is in (-pi, pi] but the task range is [0, pi]; isotropy makes |phi| equivalent
                    histograms[c].phi->Fill(std::abs(phi[i]));
                }
            }
            // once per event, also when the multiplicity is zero
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

    for (const CategoryHistograms& h : histograms)
    {
        outFile->WriteTObject(h.n.get());
        outFile->WriteTObject(h.pt.get());
        outFile->WriteTObject(h.eta.get());
        outFile->WriteTObject(h.phi.get());
    }
    outFile->WriteTObject(&generatorSettings);
    outFile->Close();

    std::cout << "Read " << nEvents << " events with " << nParticles << " final-state particles from "
              << argv[1] << "\n";
    std::cout << "Wrote " << 4 * histograms.size() << " histograms and the generator settings to " << argv[2]
              << "\n\n";

    std::cout << "Binning:\n";
    std::printf("  hN_*    multiplicity per event     %3d bins in [%g, %g]\n", nBinsN, minN, maxN);
    std::printf("  hPt_*   transverse momentum [GeV]  %3d bins in [%g, %g]\n", nBinsPt, minPt, maxPt);
    std::printf("  hEta_*  pseudorapidity             %3d bins in [%g, %g]\n", nBinsEta, minEta, maxEta);
    std::printf("  hPhi_*  azimuthal angle |phi|      %3d bins in [%g, pi]\n\n", nBinsPhi, minPhi);

    std::cout << "Particles per category:\n";
    std::printf("  %-6s %10s %10s\n", "", "in total", "per event");
    for (std::size_t c = 0; c < categories.size(); ++c)
    {
        std::printf("  %-6s %10.0f %10.2f\n", categories[c].name.c_str(), histograms[c].pt->GetEntries(),
                    histograms[c].n->GetMean());
    }

    std::cout << "\nStatistics of every histogram over all its entries; values outside the axis range\n"
                 "are kept in the underflow and overflow bins:\n";
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
