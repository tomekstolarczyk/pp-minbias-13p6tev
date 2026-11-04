// Checks the outputs of scripts/run_all.sh against the task, section by section.
// Usage, from the repository root after the pipeline: root -l -b -q tests/check_task.C
// (the event and histogram files can be given as arguments, e.g. to check another run)
// The exit code is 0 when every check passes and 1 otherwise.

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"


const std::string generateLogName = "results/generate.log";
const std::string macroFileName = "macros/draw_plots.C";
const std::string plotsDir = "plots";
const std::string compileCommandsName = "build/compile_commands.json";

// "about 100 000 events"
constexpr long long requiredEvents = 100000;

// two protons collide, so the final state of every event has total charge +2
constexpr int beamCharge = 2;

// "N: 0-400 (400 bins), pT: 0-10 GeV (100 bins), eta: -5 to 5 (100 bins), phi: 0-pi (64 bins)"
struct Binning
{
    std::string prefix;
    std::string description;
    int nBins;
    double min;
    double max;
};

const std::vector<Binning> binnings = {
    {"hN", "N: 0-400, 400 bins", 400, 0.0, 400.0},
    {"hPt", "pT: 0-10 GeV, 100 bins", 100, 0.0, 10.0},
    {"hEta", "eta: -5 to 5, 100 bins", 100, -5.0, 5.0},
    {"hPhi", "phi: 0-pi, 64 bins", 64, 0.0, TMath::Pi()},
};

// "all charged: charge != 0" and "PDG mapping: pions +-211, K +-321, p +-2212";
// "all" holds every final-state particle
struct Category
{
    std::string name;
    bool (*selects)(int pdg, int charge);
};

const std::vector<Category> categories = {
    {"all", [](int, int) { return true; }},
    {"ch", [](int, int charge) { return charge != 0; }},
    {"pip", [](int pdg, int) { return pdg == 211; }},
    {"pim", [](int pdg, int) { return pdg == -211; }},
    {"Kp", [](int pdg, int) { return pdg == 321; }},
    {"Km", [](int pdg, int) { return pdg == -321; }},
    {"p", [](int pdg, int) { return pdg == 2212; }},
    {"pbar", [](int pdg, int) { return pdg == -2212; }},
};

// "axis labels (p_{T} [GeV], \eta, \phi [rad], Multiplicity N)", as ROOT writes them
const std::vector<std::string> axisLabels = {"p_{T} [GeV]", "#eta", "#phi [rad]", "Multiplicity N"};

// prints the result of one requirement and returns it
bool report(bool passed, const std::string& requirement)
{
    std::cout << "check_task:   " << (passed ? "ok    " : "FAILED") << "  " << requirement << "\n";
    return passed;
}

// whole content of a text file, empty if it cannot be read
std::string readText(const std::string& fileName)
{
    std::ifstream file(fileName);
    std::stringstream content;
    content << file.rdbuf();
    return content.str();
}

bool contains(const std::string& text, const std::string& part)
{
    return text.find(part) != std::string::npos;
}

// the first line of a text that contains key, empty if there is none
std::string lineWith(const std::string& text, const std::string& key)
{
    const std::size_t at = text.find(key);
    if (at == std::string::npos)
    {
        return "";
    }
    const std::size_t newline = text.rfind('\n', at);
    const std::size_t begin = newline == std::string::npos ? 0 : newline + 1;
    return text.substr(begin, text.find('\n', at) - begin);
}

// section 1: the settings Pythia reports in its log and the events in the tree
bool checkGeneration(const std::string& log, long long nEvents, long long nWrongCharge)
{
    std::cout << "check_task: 1) event generation (Pythia 8)\n";
    // Pythia lists only the settings that differ from its defaults; the default beam is a proton
    const std::string beamA = lineWith(log, "Beams:idA");
    const std::string beamB = lineWith(log, "Beams:idB");
    const bool protons = (beamA.empty() || contains(beamA, "2212")) && (beamB.empty() || contains(beamB, "2212"));

    bool passed = true;
    passed &= report(protons, "pp collisions");
    passed &= report(contains(lineWith(log, "Beams:eCM"), "13600.000"), "Beams:eCM = 13600");
    passed &= report(contains(lineWith(log, "SoftQCD:"), " on "), "minimum bias (SoftQCD)");
    passed &= report(nEvents == requiredEvents, "about 100 000 events: " + std::to_string(nEvents));
    passed &= report(nWrongCharge == 0, "only final-state particles: total charge +2 in every event");
    return passed;
}

// sections 2 and 3: the histograms in analysis.root
bool checkHistograms(TFile& analysisFile, long long nEvents, const std::vector<long long>& nParticles)
{
    std::cout << "check_task: 2) and 3) histograms in analysis.root\n";
    // the file owns the histograms it returns
    auto get = [&analysisFile](const std::string& name) { return analysisFile.Get<TH1D>(name.c_str()); };

    bool allPresent = true;
    for (const Binning& binning : binnings)
    {
        for (const Category& category : categories)
        {
            allPresent = allPresent && get(binning.prefix + "_" + category.name) != nullptr;
        }
    }
    // the remaining checks read every histogram
    if (!report(allPresent, "hN_*, hPt_*, hEta_*, hPhi_* for all, ch, pip, pim, Kp, Km, p, pbar"))
    {
        return false;
    }

    bool passed = true;
    for (const Binning& binning : binnings)
    {
        bool sameBinning = true;
        for (const Category& category : categories)
        {
            const TAxis* axis = get(binning.prefix + "_" + category.name)->GetXaxis();
            sameBinning = sameBinning && axis->GetNbins() == binning.nBins && axis->GetXmin() == binning.min &&
                          axis->GetXmax() == binning.max;
        }
        passed &= report(sameBinning, binning.description);
    }

    // a multiplicity histogram gets one number from every event, also when it is zero
    bool oncePerEvent = true;
    for (const Category& category : categories)
    {
        oncePerEvent = oncePerEvent && get("hN_" + category.name)->GetEntries() == nEvents;
    }
    passed &= report(oncePerEvent, "multiplicity per event: every hN has one entry per event");

    bool phiInRange = true;
    for (const Category& category : categories)
    {
        const TH1D* phi = get("hPhi_" + category.name);
        phiInRange = phiInRange && phi->GetBinContent(0) == 0 && phi->GetBinContent(phi->GetNbinsX() + 1) == 0;
    }
    passed &= report(phiInRange, "phi: every entry in 0-pi");

    // every particle of a category is in its pT, eta and phi histograms, including under- and overflow
    bool countsMatch = true;
    for (std::size_t c = 0; c < categories.size(); ++c)
    {
        for (const std::string prefix : {"hPt", "hEta", "hPhi"})
        {
            const double entries = get(prefix + "_" + categories[c].name)->GetEntries();
            if (entries != nParticles[c])
            {
                std::cerr << prefix << "_" << categories[c].name << ": " << entries << " entries, "
                          << nParticles[c] << " particles in the tree\n";
                countsMatch = false;
            }
        }
    }
    passed &= report(countsMatch, "categories (all charged: charge != 0; PDG +-211, +-321, +-2212): "
                                  "each holds exactly its particles");
    return passed;
}

// section 4: the drawing macro and the plots it saved
bool checkDrawing()
{
    std::cout << "check_task: 4) drawing macro draw_plots.C\n";
    const std::string macro = readText(macroFileName);

    bool labels = !macro.empty();
    for (const std::string& label : axisLabels)
    {
        labels = labels && contains(macro, label);
    }

    // every histogram on its own canvas, saved as PNG and PDF
    bool singlePlots = true;
    for (const Binning& binning : binnings)
    {
        for (const Category& category : categories)
        {
            const std::string base = plotsDir + "/histograms_single/" + binning.prefix + "_" + category.name;
            singlePlots = singlePlots && std::filesystem::exists(base + ".png") &&
                          std::filesystem::exists(base + ".pdf");
        }
    }

    // the categories of one group overlaid on one canvas, told apart by a legend
    bool overlays = contains(macro, "TLegend");
    for (const std::string variable : {"multiplicity", "pt", "eta", "phi"})
    {
        for (const std::string group : {"all_ch", "pions", "kaons", "protons"})
        {
            const std::string base = plotsDir + "/histograms_overlay/" + variable + "_" + group;
            overlays = overlays && std::filesystem::exists(base + ".png") && std::filesystem::exists(base + ".pdf");
        }
    }

    bool passed = true;
    passed &= report(contains(macro, "analysis.root"), "opens analysis.root");
    passed &= report(labels, "axis labels p_{T} [GeV], #eta, #phi [rad], Multiplicity N");
    passed &= report(singlePlots, "every histogram on its own canvas, exported as PNG and PDF");
    passed &= report(overlays, "overlays with legends, exported as PNG and PDF");
    return passed;
}

// technical requirements: the C++ standard of our build and the versions of ROOT and Pythia
bool checkTechnical(const std::string& log)
{
    std::cout << "check_task: technical requirements\n";
    const std::string commands = readText(compileCommandsName);
    const bool cpp17 = contains(commands, "-std=c++17") && !contains(commands, "-std=c++2") &&
                       !contains(commands, "-std=gnu++");
    const std::string pythiaLine = lineWith(log, "This is PYTHIA version 8.");

    bool passed = true;
    passed &= report(cpp17, "C++17: generate and analysis compiled with -std=c++17");
    passed &= report(gROOT->GetVersionInt() / 10000 == 6, std::string("ROOT 6: ") + gROOT->GetVersion());
    passed &= report(!pythiaLine.empty(), "Pythia 8: the log says it was PYTHIA 8");
    return passed;
}

int check_task(const char* eventsFileName = "data/events.root",
               const char* analysisFileName = "results/analysis.root")
{
    std::unique_ptr<TFile> eventsFile(TFile::Open(eventsFileName, "READ"));
    std::unique_ptr<TFile> analysisFile(TFile::Open(analysisFileName, "READ"));
    if (!eventsFile || !analysisFile)
    {
        std::cerr << "Could not open " << eventsFileName << " or " << analysisFileName << "\n";
        return 1;
    }
    TTree* tree = eventsFile->Get<TTree>("events");
    if (!tree)
    {
        std::cerr << "No tree \"events\" in " << eventsFileName << "\n";
        return 1;
    }

    // read the tree once
    long long nEvents = 0;
    long long nWrongCharge = 0;
    std::vector<long long> nParticles(categories.size(), 0);
    {
        // the reader uses the tree in its destructor, so this block ends before the file is closed
        TTreeReader reader(tree);
        TTreeReaderArray<int> pdg(reader, "pdg");
        TTreeReaderArray<int> charge(reader, "charge");
        while (reader.Next())
        {
            ++nEvents;
            int totalCharge = 0;
            for (std::size_t i = 0; i < pdg.GetSize(); ++i)
            {
                totalCharge += charge[i];
                for (std::size_t c = 0; c < categories.size(); ++c)
                {
                    if (categories[c].selects(pdg[i], charge[i]))
                    {
                        ++nParticles[c];
                    }
                }
            }
            if (totalCharge != beamCharge)
            {
                ++nWrongCharge;
            }
        }
        if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd)
        {
            std::cerr << "Error while reading " << eventsFileName << " at event " << nEvents << "\n";
            return 1;
        }
    }

    const std::string log = readText(generateLogName);
    bool passed = true;
    passed &= checkTechnical(log);
    passed &= checkGeneration(log, nEvents, nWrongCharge);
    passed &= checkHistograms(*analysisFile, nEvents, nParticles);
    passed &= checkDrawing();

    std::cout << "check_task: " << (passed ? "PASSED" : "FAILED") << "\n";
    return passed ? 0 : 1;
}
