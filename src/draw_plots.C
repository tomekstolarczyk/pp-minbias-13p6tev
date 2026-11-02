// Draw the histograms from analysis.root and save every plot as PNG and PDF
//
// Usage: root -l -b -q 'src/draw_plots.C("results/analysis.root", "plots")'
// The exit code is 0 when every plot was saved and 1 otherwise.

#include <algorithm>
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
#include "TLegend.h"
#include "TStyle.h"

struct Curve
{
    std::string category;
    std::string legend;
    Color_t color;
    Style_t lineStyle;
    Style_t markerStyle;
};

struct Variable
{
    std::string prefix;
    std::string fileTag;
    std::string title;
    std::string xTitle;
    std::string yTitle;
    bool logY;
    bool points;
};

// draw all curves of one variable on one canvas and save it as PNG and PDF
bool drawPlot(TFile& file, const Variable& variable, const std::vector<Curve>& curves, const std::string& groupTitle, const std::string& groupTag, const std::string& outDir)
{
    std::vector<TH1D*> histograms;
    double maximum = 0.0;
    for (const Curve& curve : curves)
    {
        const std::string name = variable.prefix + "_" + curve.category;
        TH1D* histogram = file.Get<TH1D>(name.c_str());
        if (!histogram)
        {
            std::cerr << "No histogram " << name << " in " << file.GetName() << "\n";
            return false;
        }
        histogram->SetLineColor(curve.color);
        histogram->SetLineStyle(curve.lineStyle);
        histogram->SetLineWidth(2);
        histogram->SetMarkerColor(curve.color);
        histogram->SetMarkerStyle(curve.markerStyle);
        histogram->SetMarkerSize(0.6);
        maximum = std::max(maximum, histogram->GetMaximum());
        histograms.push_back(histogram);
    }

    // the first histogram sets the frame: title, axis titles and y range for all curves
    TH1D* frame = histograms.front();
    frame->SetTitle((variable.title + ", " + groupTitle + ";" + variable.xTitle + ";" + variable.yTitle).c_str());
    if (variable.logY)
    {
        // counts are whole numbers, so 0.5 shows every non-empty bin on a log scale
        frame->SetMinimum(0.5);
        frame->SetMaximum(maximum * 50.0);
    }
    else
    {
        // start at zero
        frame->SetMinimum(0.0);
        frame->SetMaximum(maximum * 1.5);
    }
    frame->GetYaxis()->SetTitleOffset(1.4);

    // particle and antiparticle share a row, so six entries take three rows
    const int nColumns = curves.size() > 3 ? 2 : 1;
    const int nRows = (static_cast<int>(curves.size()) + nColumns - 1) / nColumns;
    TLegend legend(0.62, 0.86 - 0.06 * nRows, 0.88, 0.86);
    legend.SetNColumns(nColumns);
    TCanvas canvas(("c_" + variable.fileTag + "_" + groupTag).c_str(), groupTitle.c_str(), 800, 600);
    canvas.SetLeftMargin(0.13);
    canvas.SetLogy(variable.logY);

    const std::string drawOption = variable.points ? "P" : "HIST";
    for (std::size_t i = 0; i < histograms.size(); ++i)
    {
        histograms[i]->Draw(i == 0 ? drawOption.c_str() : (drawOption + " SAME").c_str());
        legend.AddEntry(histograms[i], curves[i].legend.c_str(), variable.points ? "p" : "l");
    }
    legend.Draw();

    const std::string baseName = outDir + "/" + variable.fileTag + "_" + groupTag;
    canvas.SaveAs((baseName + ".png").c_str());
    canvas.SaveAs((baseName + ".pdf").c_str());
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

    // plots with several histograms show the legend instead of the statistics box
    gStyle->SetOptStat(0);
    // axis labels get at most three digits; larger numbers share a common factor such as x10^3
    TGaxis::SetMaxDigits(3);

    const std::vector<Variable> variables = {
        {"hN", "multiplicity", "Multiplicity", "Multiplicity N", "Number of events", true, true},
        {"hPt", "pt", "Transverse momentum", "p_{T} [GeV]", "Number of particles", true, false},
        {"hEta", "eta", "Pseudorapidity", "#eta", "Number of particles", false, false},
        {"hPhi", "phi", "Azimuthal angle", "#phi [rad]", "Number of particles", false, false},
    };

    const std::vector<Curve> inclusive = {
        {"all", "all final-state", kBlack, kSolid, kFullCircle},
        {"ch", "charged", kMagenta + 2, kSolid, kFullCircle},
    };
    // colour tells the particle species; solid lines and full markers mean positive charge
    const std::vector<Curve> identified = {
        {"pip", "#pi^{+}", kRed + 1, kSolid, kFullCircle},
        {"pim", "#pi^{-}", kRed + 1, kDashed, kOpenCircle},
        {"Kp", "K^{+}", kBlue + 1, kSolid, kFullCircle},
        {"Km", "K^{-}", kBlue + 1, kDashed, kOpenCircle},
        {"p", "p", kGreen + 2, kSolid, kFullCircle},
        {"pbar", "#bar{p}", kGreen + 2, kDashed, kOpenCircle},
    };

    for (const Variable& variable : variables)
    {
        if (!drawPlot(*file, variable, inclusive, "all and charged particles", "all_ch", outDir) ||
            !drawPlot(*file, variable, identified, "identified hadrons", "identified", outDir))
        {
            return 1;
        }
    }
    return 0;
}
