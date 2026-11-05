// Compares the charged-particle density dN_ch/deta with the ALICE measurement in pp collisions at
// sqrt(s) = 13.6 TeV (ALICE, arXiv:2602.10658). This is an extra, not part of the task.
// Usage: root -l -b -q 'macros/compare_alice.C("data/events.root", "data/events_alice_primaries.root", "plots/alice_comparison")'
// The exit code is 0 when both files were read and the plots saved, 1 otherwise.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <system_error>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

// ALICE definitions (arXiv:2602.10658)
// INEL>0: at least one charged particle in |eta| < 1
constexpr double inelGt0EtaMax = 1.0;
// the density is averaged over |eta| < 0.5, a window of width 1, so it equals the mean count there
constexpr double midEtaMax = 0.5;
// multiplicity classes use the forward detector FT0: FT0A 3.5 < eta < 4.9, FT0C -3.3 < eta < -2.1
constexpr double ft0aMin = 3.5;
constexpr double ft0aMax = 4.9;
constexpr double ft0cMin = -3.3;
constexpr double ft0cMax = -2.1;

// measured <dN_ch/deta> in |eta| < 0.5 for INEL>0 events, all uncertainties combined
constexpr double aliceDensity = 7.10;
constexpr double aliceDensityError = 0.18;

// Table 2 of the paper: multiplicity classes as INEL>0 percentiles (highest multiplicity first),
// measured <dN_ch/deta> and the PYTHIA 8.311 Monash value computed by ALICE
struct MultiplicityClass
{
    std::string label; // percentile range used by ALICE
    double fromPercent;
    double toPercent;
    double data;
    double dataError;
    double alicePythia;
};

const std::vector<MultiplicityClass> multiplicityClasses = {
    {"0-1%", 0.0, 0.96, 21.78, 0.38, 25.50},     {"1-5%", 0.96, 4.80, 18.48, 0.25, 20.71},
    {"5-10%", 4.80, 9.60, 15.76, 0.22, 17.19},   {"10-15%", 9.60, 14.42, 13.89, 0.19, 14.71},
    {"15-20%", 14.42, 19.25, 12.50, 0.17, 12.62}, {"20-30%", 19.25, 28.92, 10.86, 0.15, 10.02},
    {"30-40%", 28.92, 38.62, 9.09, 0.13, 7.51},   {"40-50%", 38.62, 48.38, 7.63, 0.11, 5.60},
    {"50-70%", 48.38, 68.09, 5.87, 0.09, 3.72},   {"70-100%", 68.09, 100.0, 3.69, 0.06, 2.47},
};

// eta binning of the density plot, as in the task
constexpr int nBinsEta = 100;
constexpr double minEta = -5.0;
constexpr double maxEta = 5.0;

constexpr double textSize = 0.035;

// mean of a count and the statistical error of that mean
struct Mean
{
    double value;
    double error;
};

Mean meanOf(double sum, double sumOfSquares, double n)
{
    const double mean = sum / n;
    const double variance = sumOfSquares / n - mean * mean;
    return {mean, std::sqrt(variance / n)};
}

// what the multiplicity classes need from one INEL>0 event
struct InelGt0Event
{
    int nForward; // charged particles in the FT0 acceptance, used to sort the events
    int nMid;     // charged particles in |eta| < 0.5
};

// what one run gives: event counts, charged particles in |eta| < 0.5 summed over all events and over
// INEL>0 events, eta distributions of charged particles, and the INEL>0 events for the classes
struct Run
{
    long long nEvents = 0;
    long long nInelGt0 = 0;
    double sumMid = 0.0;
    double sumMidSquares = 0.0;
    double sumMidInelGt0 = 0.0;
    double sumMidInelGt0Squares = 0.0;
    std::unique_ptr<TH1D> etaAll;
    std::unique_ptr<TH1D> etaInelGt0;
    std::vector<InelGt0Event> inelGt0Events;
};

std::unique_ptr<TH1D> makeEtaHistogram(const std::string& name)
{
    auto histogram = std::make_unique<TH1D>(name.c_str(), "", nBinsEta, minEta, maxEta);
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    return histogram;
}

bool inForwardDetector(double eta)
{
    return (eta > ft0aMin && eta < ft0aMax) || (eta > ft0cMin && eta < ft0cMax);
}

bool readRun(const std::string& fileName, const std::string& tag, Run& run)
{
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    if (!file)
    {
        std::cerr << "Could not open " << fileName << "\n";
        return false;
    }
    TTree* tree = file->Get<TTree>("events");
    if (!tree)
    {
        std::cerr << "No tree \"events\" in " << fileName << "\n";
        return false;
    }
    run.etaAll = makeEtaHistogram("etaAll_" + tag);
    run.etaInelGt0 = makeEtaHistogram("etaInelGt0_" + tag);

    // the reader uses the tree in its destructor, so this block ends before the file is closed
    {
        TTreeReader reader(tree);
        TTreeReaderArray<int> charge(reader, "charge");
        TTreeReaderArray<double> eta(reader, "eta");
        std::vector<double> chargedEta;
        while (reader.Next())
        {
            ++run.nEvents;
            chargedEta.clear();
            int nMid = 0;
            int nForward = 0;
            bool inelGt0 = false;
            for (std::size_t i = 0; i < charge.GetSize(); ++i)
            {
                if (charge[i] == 0)
                {
                    continue;
                }
                chargedEta.push_back(eta[i]);
                inelGt0 = inelGt0 || std::abs(eta[i]) < inelGt0EtaMax;
                nMid += std::abs(eta[i]) < midEtaMax ? 1 : 0;
                nForward += inForwardDetector(eta[i]) ? 1 : 0;
            }
            run.sumMid += nMid;
            run.sumMidSquares += static_cast<double>(nMid) * nMid;
            for (double value : chargedEta)
            {
                run.etaAll->Fill(value);
            }
            if (inelGt0)
            {
                ++run.nInelGt0;
                run.sumMidInelGt0 += nMid;
                run.sumMidInelGt0Squares += static_cast<double>(nMid) * nMid;
                for (double value : chargedEta)
                {
                    run.etaInelGt0->Fill(value);
                }
                run.inelGt0Events.push_back({nForward, nMid});
            }
        }
        if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd)
        {
            std::cerr << "Error while reading " << fileName << "\n";
            return false;
        }
    }

    // particles per event and per unit of eta; Sumw2 scales the bin errors with the contents
    run.etaAll->Scale(1.0 / (run.nEvents * run.etaAll->GetBinWidth(1)));
    run.etaInelGt0->Scale(1.0 / (run.nInelGt0 * run.etaInelGt0->GetBinWidth(1)));
    return true;
}

// sorts the INEL>0 events by the number of charged particles in the FT0 acceptance, highest first, so
// that a percentile range is a range of indices; events with the same forward count keep their
// generation order, so the split between them is arbitrary
void sortByForwardActivity(std::vector<InelGt0Event>& events)
{
    std::stable_sort(events.begin(), events.end(),
                     [](const InelGt0Event& a, const InelGt0Event& b) { return a.nForward > b.nForward; });
}

// index of the sorted event at a percentile; a class is the index range [from, to)
std::size_t percentileIndex(const std::vector<InelGt0Event>& events, double percent)
{
    return static_cast<std::size_t>(std::llround(events.size() * percent / 100.0));
}

// <dN_ch/deta> in |eta| < 0.5 of the sorted events between two percentiles
Mean densityInRange(const std::vector<InelGt0Event>& events, double fromPercent, double toPercent)
{
    const std::size_t first = percentileIndex(events, fromPercent);
    const std::size_t last = percentileIndex(events, toPercent);
    double sum = 0.0;
    double sumOfSquares = 0.0;
    for (std::size_t i = first; i < last; ++i)
    {
        sum += events[i].nMid;
        sumOfSquares += static_cast<double>(events[i].nMid) * events[i].nMid;
    }
    return meanOf(sum, sumOfSquares, static_cast<double>(last - first));
}

// text at a position given as a fraction of the canvas
void drawLabel(double x, double y, const std::string& text)
{
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(textSize);
    latex.DrawLatex(x, y, text.c_str());
}

void save(TCanvas& canvas, const std::string& baseName)
{
    canvas.SaveAs((baseName + ".png").c_str());
    canvas.SaveAs((baseName + ".pdf").c_str());
}

// Fig. 1 of the paper: dN_ch/deta (per event) for the three steps of the anatomy, with the ALICE value
void drawDensity(const Run& base, const Run& primaries, const std::string& outDir)
{
    TH1D& baseAll = *base.etaAll;
    TH1D& baseInelGt0 = *base.etaInelGt0;
    TH1D& primariesInelGt0 = *primaries.etaInelGt0;

    // dN_ch/deta means per event, as in the paper
    baseAll.SetTitle("Charged-particle density;#eta;dN_{ch}/d#eta");
    // the range of Fig. 1, where ALICE measures; the full eta range is in the task plots
    baseAll.GetXaxis()->SetRangeUser(-1.5, 1.5);
    baseAll.SetMinimum(5.0);
    baseAll.SetMaximum(9.5);
    baseAll.SetStats(false);
    baseAll.GetYaxis()->SetTitleOffset(1.2);
    // the curve with the ALICE definitions stands out, the two intermediate steps are thinner
    baseAll.SetLineColor(kGray + 1);
    baseAll.SetLineStyle(kDotted);
    baseAll.SetLineWidth(2);
    baseInelGt0.SetLineColor(kAzure + 2);
    baseInelGt0.SetLineStyle(kDashed);
    baseInelGt0.SetLineWidth(2);
    primariesInelGt0.SetLineColor(kRed + 1);
    primariesInelGt0.SetLineWidth(3);

    // the ALICE value is an average over |eta| < 0.5, drawn at eta = 0 with its total uncertainty
    TGraphErrors alice(1);
    alice.SetPoint(0, 0.0, aliceDensity);
    alice.SetPointError(0, 0.0, aliceDensityError);
    alice.SetMarkerStyle(kFullCircle);
    alice.SetMarkerSize(1.2);
    alice.SetLineWidth(2);

    TLegend legend(0.15, 0.13, 0.85, 0.37);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(textSize);
    legend.AddEntry(&alice, "ALICE data, INEL>0, mean over #left|#eta#right| < 0.5", "pe");
    legend.AddEntry(&primariesInelGt0, "PYTHIA, as ALICE: INEL>0, primary particles", "l");
    legend.AddEntry(&baseInelGt0, "PYTHIA, INEL>0, with weak-decay products", "l");
    legend.AddEntry(&baseAll, "PYTHIA, all inelastic, with weak-decay products", "l");

    TCanvas canvas("c_density", "", 800, 600);
    canvas.SetLeftMargin(0.13);
    canvas.SetTicks(1, 1);
    baseAll.Draw("HIST");
    baseInelGt0.Draw("HIST SAME");
    primariesInelGt0.Draw("HIST SAME");
    alice.Draw("P SAME");
    legend.Draw();
    drawLabel(0.17, 0.84, "PYTHIA 8 (Monash 2013), pp, #sqrt{s} = 13.6 TeV");
    save(canvas, outDir + "/alice_fig1_dndeta");
}

int compare_alice(const char* baseFileName = "data/events.root",
                  const char* primariesFileName = "data/events_alice_primaries.root",
                  const char* outDir = "plots/alice_comparison")
{
    Run base;
    Run primaries;
    if (!readRun(baseFileName, "base", base) || !readRun(primariesFileName, "primaries", primaries))
    {
        return 1;
    }
    std::error_code error;
    std::filesystem::create_directories(outDir, error);
    if (error)
    {
        std::cerr << "Could not create " << outDir << ": " << error.message() << "\n";
        return 1;
    }

    // 1. anatomy of the comparison
    const Mean baseInel = meanOf(base.sumMid, base.sumMidSquares, base.nEvents);
    const Mean baseInelGt0 = meanOf(base.sumMidInelGt0, base.sumMidInelGt0Squares, base.nInelGt0);
    const Mean primariesInel = meanOf(primaries.sumMid, primaries.sumMidSquares, primaries.nEvents);
    const Mean primariesInelGt0 =
        meanOf(primaries.sumMidInelGt0, primaries.sumMidInelGt0Squares, primaries.nInelGt0);

    std::printf("Comparison with ALICE, pp at 13.6 TeV (arXiv:2602.10658)\n\n");
    std::printf("INEL>0 events: task run %lld of %lld (%.2f%%), primaries run %lld of %lld (%.2f%%)\n\n",
                base.nInelGt0, base.nEvents, 100.0 * base.nInelGt0 / base.nEvents, primaries.nInelGt0,
                primaries.nEvents, 100.0 * primaries.nInelGt0 / primaries.nEvents);
    std::printf("<dN_ch/deta> in |eta| < 0.5, step by step towards the ALICE definitions:\n");
    std::printf("  %-52s %6.3f +- %.3f\n", "task run: INEL, default decays", baseInel.value, baseInel.error);
    std::printf("  %-52s %6.3f +- %.3f\n", "task run: INEL>0, default decays", baseInelGt0.value,
                baseInelGt0.error);
    std::printf("  %-52s %6.3f +- %.3f\n", "primaries run: INEL, ALICE primary particles", primariesInel.value,
                primariesInel.error);
    std::printf("  %-52s %6.3f +- %.3f\n", "primaries run: INEL>0, ALICE primary particles",
                primariesInelGt0.value, primariesInelGt0.error);
    std::printf("  %-52s %6.3f +- %.3f\n", "ALICE: INEL>0, primary particles", aliceDensity,
                aliceDensityError);
    std::printf("  matched definitions / ALICE = %.3f, difference = %.1f sigma (ALICE uncertainty)\n\n",
                primariesInelGt0.value / aliceDensity, (primariesInelGt0.value - aliceDensity) / aliceDensityError);

    // 3. multiplicity classes: Table 2 with our PYTHIA
    sortByForwardActivity(primaries.inelGt0Events);
    std::printf("<dN_ch/deta> in |eta| < 0.5 per FT0 multiplicity class (INEL>0, ALICE primary particles):\n");
    std::printf("  %-8s %15s %14s %16s %12s\n", "class", "ALICE data", "ALICE PYTHIA", "our PYTHIA", "ours / data");
    for (std::size_t i = 0; i < multiplicityClasses.size(); ++i)
    {
        const MultiplicityClass& c = multiplicityClasses[i];
        const Mean ours = densityInRange(primaries.inelGt0Events, c.fromPercent, c.toPercent);
        std::printf("  %-8s %8.2f +- %.2f %14.2f %9.2f +- %.2f %12.3f\n", c.label.c_str(), c.data, c.dataError,
                    c.alicePythia, ours.value, ours.error, ours.value / c.data);
    }

    // 2. plot
    gStyle->SetOptStat(0);
    drawDensity(base, primaries, outDir);
    return 0;
}
