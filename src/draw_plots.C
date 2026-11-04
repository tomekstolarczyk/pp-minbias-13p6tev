// Draws the histograms from analysis.root and saves every plot as PNG and PDF:
// each histogram on its own canvas (32 plots) and, for every variable, the categories
// of one group overlaid with a legend (8 plots).
// Usage: root -l -b -q 'src/draw_plots.C("results/analysis.root", "plots")'
// The exit code is 0 when every plot was saved and 1 otherwise.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <system_error>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

// one of the four histogrammed quantities and how to draw it
struct Variable
{
    std::string prefix;  // histogram name prefix, e.g. hPt
    std::string fileTag; // used in the file names of the overlay plots, e.g. pt
    std::string title;
    std::string xTitle;
    std::string yTitle;
    bool logY;   // logarithmic y axis, for counts spanning several orders of magnitude
    bool points; // markers with error bars instead of a line
};

// one category drawn as a curve
struct Curve
{
    std::string category; // histogram name suffix, e.g. pip
    std::string legend;   // text in the legend of an overlay plot
    std::string title;    // used in the title of the single plot
    Color_t color;
    Style_t lineStyle;
    Style_t markerStyle;
};

// categories that are compared on one overlay plot
struct Group
{
    std::string title;
    std::string tag; // used in the file names, e.g. identified
    std::vector<Curve> curves;
};

// text size of the labels and legends
constexpr double textSize = 0.035;

// reads one histogram from the file and sets its colour, line and marker; the file owns it
TH1D* getStyledHistogram(TFile& file, const Variable& variable, const Curve& curve)
{
    const std::string name = variable.prefix + "_" + curve.category;
    TH1D* histogram = file.Get<TH1D>(name.c_str());
    if (!histogram)
    {
        std::cerr << "No histogram " << name << " in " << file.GetName() << "\n";
        return nullptr;
    }

    histogram->SetLineColor(curve.color);
    histogram->SetLineStyle(variable.points ? kSolid : curve.lineStyle);
    histogram->SetLineWidth(2);
    histogram->SetMarkerColor(curve.color);
    histogram->SetMarkerStyle(curve.markerStyle);
    histogram->SetMarkerSize(0.6);
    return histogram;
}

// content of the highest bin; GetMaximum() would return a value set earlier with SetMaximum()
double highestBin(const TH1D& histogram)
{
    return histogram.GetBinContent(histogram.GetMaximumBin());
}

// the first histogram drawn on a canvas sets the frame: title, axis titles and y range
void setFrame(TH1D& frame, const Variable& variable, const std::string& title, double maximum)
{
    frame.SetTitle((title + ";" + variable.xTitle + ";" + variable.yTitle).c_str());
    frame.GetYaxis()->SetTitleOffset(1.4);

    if (variable.logY)
    {
        // counts are whole numbers, so 0.5 shows every non-empty bin on a log scale
        frame.SetMinimum(0.5);
        frame.SetMaximum(maximum * 50.0);
    }
    else
    {
        // without this ROOT starts a linear axis just below the lowest bin
        frame.SetMinimum(0.0);
        frame.SetMaximum(maximum * 1.5);
    }
}

void prepareCanvas(TCanvas& canvas, const Variable& variable)
{
    canvas.SetLeftMargin(0.13);
    canvas.SetLogy(variable.logY);
    // tick marks on all four sides
    canvas.SetTicks(1, 1);
}

// tells where the data come from
void drawInfo(const std::string& nEvents)
{
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(textSize);
    latex.DrawLatex(0.17, 0.84, "PYTHIA 8 (Monash 2013)");
    latex.DrawLatex(0.17, 0.79, "pp, #sqrt{s} = 13.6 TeV");
    latex.DrawLatex(0.17, 0.74, ("minimum bias, " + nEvents + " events").c_str());
}

void save(TCanvas& canvas, const std::string& baseName)
{
    canvas.SaveAs((baseName + ".png").c_str());
    canvas.SaveAs((baseName + ".pdf").c_str());
}

// one histogram on its own canvas; the statistics box shows what lies outside the range
bool drawSingle(TFile& file, const Variable& variable, const Curve& curve, const std::string& nEvents,
                const std::string& outDir)
{
    TH1D* histogram = getStyledHistogram(file, variable, curve);
    if (!histogram)
    {
        return false;
    }
    histogram->SetStats(true);
    histogram->SetLineStyle(kSolid);
    setFrame(*histogram, variable, variable.title + ", " + curve.title, highestBin(*histogram));

    TCanvas canvas((std::string("c_") + histogram->GetName()).c_str(), "", 800, 600);
    prepareCanvas(canvas, variable);
    histogram->Draw(variable.points ? "P" : "HIST");
    drawInfo(nEvents);
    save(canvas, outDir + "/" + histogram->GetName());
    return true;
}

// all categories of one group on one canvas, told apart by the legend
bool drawOverlay(TFile& file, const Variable& variable, const Group& group, const std::string& nEvents,
                 const std::string& outDir)
{
    std::vector<TH1D*> histograms;
    double maximum = 0.0;
    for (const Curve& curve : group.curves)
    {
        TH1D* histogram = getStyledHistogram(file, variable, curve);
        if (!histogram)
        {
            return false;
        }
        // the legend takes the place of the statistics box
        histogram->SetStats(false);
        maximum = std::max(maximum, highestBin(*histogram));
        histograms.push_back(histogram);
    }
    setFrame(*histograms.front(), variable, variable.title + ", " + group.title, maximum);

    // particle and antiparticle share a row, so six entries take three rows
    const int nColumns = group.curves.size() > 3 ? 2 : 1;
    const int nRows = (static_cast<int>(group.curves.size()) + nColumns - 1) / nColumns;
    // the canvas still lists the legend when it is destroyed, so the legend is declared first
    TLegend legend(0.62, 0.86 - 0.06 * nRows, 0.88, 0.86);
    legend.SetNColumns(nColumns);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(textSize);
    TCanvas canvas(("c_" + variable.fileTag + "_" + group.tag).c_str(), "", 800, 600);
    prepareCanvas(canvas, variable);

    const std::string option = variable.points ? "P" : "HIST";
    for (std::size_t i = 0; i < histograms.size(); ++i)
    {
        // SAME draws into the frame of the first histogram
        histograms[i]->Draw(i == 0 ? option.c_str() : (option + " SAME").c_str());
        legend.AddEntry(histograms[i], group.curves[i].legend.c_str(), variable.points ? "p" : "l");
    }
    legend.Draw();
    drawInfo(nEvents);
    save(canvas, outDir + "/" + variable.fileTag + "_" + group.tag);
    return true;
}

int draw_plots(const char* inFileName = "results/analysis.root", const char* outDir = "plots")
{
    std::unique_ptr<TFile> file(TFile::Open(inFileName, "READ"));
    if (!file)
    {
        std::cerr << "Could not open " << inFileName << "\n";
        return 1;
    }

    std::error_code error;
    std::filesystem::create_directories(outDir, error);
    if (error)
    {
        std::cerr << "Could not create " << outDir << ": " << error.message() << "\n";
        return 1;
    }

    // statistics box: entries, underflow and overflow
    gStyle->SetOptStat("eou");
    gStyle->SetStatX(0.88);
    gStyle->SetStatY(0.88);
    gStyle->SetStatFormat("8.0f");
    gStyle->SetTitleFontSize(0.045);
    // axis labels get at most three digits; larger numbers share a common factor such as x10^3
    TGaxis::SetMaxDigits(3);

    // every multiplicity histogram has exactly one entry per event
    const TH1D* nCharged = file->Get<TH1D>("hN_ch");
    if (!nCharged)
    {
        std::cerr << "No histogram hN_ch in " << inFileName << "\n";
        return 1;
    }
    const std::string nEvents = std::to_string(std::llround(nCharged->GetEntries()));

    const std::vector<Variable> variables = {
        {"hN", "multiplicity", "Multiplicity", "Multiplicity N", "Events", true, true},
        {"hPt", "pt", "Transverse momentum", "p_{T} [GeV]", "Particles / 0.1 GeV", true, false},
        {"hEta", "eta", "Pseudorapidity", "#eta", "Particles / 0.1", false, false},
        {"hPhi", "phi", "Azimuthal angle", "#phi [rad]", "Particles / (#pi/64 rad)", false, false},
    };

    // colour tells the particle species; solid lines and full markers mean positive charge
    const std::vector<Group> groups = {
        {"all and charged particles", "all_ch",
         {
             {"all", "all final-state", "all final-state particles", kBlack, kSolid, kFullCircle},
             {"ch", "charged", "charged particles", kMagenta + 2, kSolid, kFullCircle},
         }},
        {"identified hadrons", "identified",
         {
             {"pip", "#pi^{+}", "#pi^{+}", kRed + 1, kSolid, kFullCircle},
             {"pim", "#pi^{-}", "#pi^{-}", kRed + 1, kDashed, kOpenCircle},
             {"Kp", "K^{+}", "K^{+}", kBlue + 1, kSolid, kFullCircle},
             {"Km", "K^{-}", "K^{-}", kBlue + 1, kDashed, kOpenCircle},
             {"p", "p", "p", kGreen + 2, kSolid, kFullCircle},
             {"pbar", "#bar{p}", "#bar{p}", kGreen + 2, kDashed, kOpenCircle},
         }},
    };

    for (const Variable& variable : variables)
    {
        for (const Group& group : groups)
        {
            if (!drawOverlay(*file, variable, group, nEvents, outDir))
            {
                return 1;
            }
            for (const Curve& curve : group.curves)
            {
                if (!drawSingle(*file, variable, curve, nEvents, outDir))
                {
                    return 1;
                }
            }
        }
    }
    return 0;
}
