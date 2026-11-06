// Draws the histograms from analysis.root and saves every plot as PNG and PDF
// Usage: root -l -b -q 'macros/draw_plots.C("results/analysis.root", "plots")'

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
#include "TLine.h"
#include "TStyle.h"

struct Variable
{
    std::string prefix;
    std::string fileTag;
    std::string title;
    std::string xTitle;
    std::string yTitle;
    bool logY;
    bool points;
    bool ratio;
};

struct Curve
{
    std::string category;
    std::string legend;
    std::string title;
    Color_t color;
    Style_t lineStyle;
    Style_t markerStyle;
};

struct Group
{
    std::string title;
    std::string tag;
    std::vector<Curve> curves;
    std::string ratioTitle;
    double ratioMin;
    double ratioMax;
    bool subset;
};

constexpr double textSize = 0.035;
constexpr double ratioPadHeight = 0.3;

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

// GetMaximum() would return a value set earlier with SetMaximum()
double highestBin(const TH1D& histogram)
{
    return histogram.GetBinContent(histogram.GetMaximumBin());
}

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

void prepareCanvas(TVirtualPad& canvas, const Variable& variable)
{
    canvas.SetLeftMargin(0.13);
    canvas.SetLogy(variable.logY);
    canvas.SetTicks(1, 1);
}

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

std::unique_ptr<TH1D> makeRatio(const TH1D& first, const TH1D& second, const Variable& variable,
                                const Group& group)
{
    std::unique_ptr<TH1D> ratio(static_cast<TH1D*>(second.Clone()));
    ratio->SetDirectory(nullptr);
    // binomial errors when one category is part of the other (charged of all)
    ratio->Divide(&second, &first, 1.0, 1.0, group.subset ? "B" : "");
    ratio->SetTitle(("; " + variable.xTitle + ";" + group.ratioTitle).c_str());
    ratio->SetMinimum(group.ratioMin);
    ratio->SetMaximum(group.ratioMax);
    ratio->SetLineColor(kBlack);
    ratio->SetLineStyle(kSolid);
    ratio->SetMarkerColor(kBlack);
    ratio->SetMarkerStyle(kFullCircle);
    ratio->SetMarkerSize(0.5);

    const double scale = (1.0 - ratioPadHeight) / ratioPadHeight;
    TAxis* x = ratio->GetXaxis();
    TAxis* y = ratio->GetYaxis();
    x->SetLabelSize(textSize * scale);
    x->SetTitleSize(textSize * scale);
    y->SetLabelSize(textSize * scale);
    y->SetTitleSize(textSize * scale);
    y->SetTitleOffset(1.4 / scale);
    y->SetNdivisions(505);
    return ratio;
}

bool drawOverlay(TFile& file, const Variable& variable, const Group& group, const std::string& nEvents,
                 const std::string& outDir)
{
    std::vector<std::unique_ptr<TH1D>> histograms;
    double maximum = 0.0;
    int lastFilledBin = 1;
    for (const Curve& curve : group.curves)
    {
        const TH1D* original = getStyledHistogram(file, variable, curve);
        if (!original)
        {
            return false;
        }
        std::unique_ptr<TH1D> histogram(static_cast<TH1D*>(original->Clone()));
        histogram->SetDirectory(nullptr);
        histogram->SetStats(false);
        maximum = std::max(maximum, highestBin(*histogram));
        lastFilledBin = std::max(lastFilledBin, histogram->FindLastBinAbove(0.0));
        histograms.push_back(std::move(histogram));
    }
    TH1D& frame = *histograms.front();
    setFrame(frame, variable, variable.title + ", " + group.title, maximum);
    // the single plots show the full range from the task; an overlay ends at the last filled bin
    frame.GetXaxis()->SetRange(1, lastFilledBin);

    // objects drawn on the canvas are declared before it
    std::unique_ptr<TH1D> ratio;
    if (variable.ratio)
    {
        ratio = makeRatio(*histograms[0], *histograms[1], variable, group);
        ratio->GetXaxis()->SetRange(1, lastFilledBin);
    }
    TLine one(frame.GetXaxis()->GetBinLowEdge(1), 1.0, frame.GetXaxis()->GetBinUpEdge(lastFilledBin), 1.0);
    one.SetLineStyle(kDashed);
    one.SetLineColor(kGray + 2);
    TLegend legend(0.62, 0.86 - 0.06 * group.curves.size(), 0.88, 0.86);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(textSize);
    TCanvas canvas(("c_" + variable.fileTag + "_" + group.tag).c_str(), "", 800, variable.ratio ? 800 : 600);

    // with a ratio the canvas is split into two pads owned by the canvas
    TVirtualPad* plotPad = &canvas;
    if (variable.ratio)
    {
        canvas.Divide(1, 2);
        plotPad = canvas.cd(1);
        plotPad->SetPad(0.0, ratioPadHeight, 1.0, 1.0);
        plotPad->SetBottomMargin(0.02);
        frame.GetXaxis()->SetLabelSize(0.0);
        frame.GetXaxis()->SetTitleSize(0.0);

        TVirtualPad* ratioPad = canvas.cd(2);
        ratioPad->SetPad(0.0, 0.0, 1.0, ratioPadHeight);
        ratioPad->SetTopMargin(0.03);
        ratioPad->SetBottomMargin(0.35);
        ratioPad->SetLeftMargin(0.13);
        ratioPad->SetTicks(1, 1);
        ratioPad->SetGridy();
        ratio->Draw("E");
        if (group.ratioMin < 1.0 && group.ratioMax > 1.0)
        {
            one.Draw();
        }
    }

    plotPad->cd();
    prepareCanvas(*plotPad, variable);
    const std::string option = variable.points ? "P" : "HIST";
    for (std::size_t i = 0; i < histograms.size(); ++i)
    {
        histograms[i]->Draw(i == 0 ? option.c_str() : (option + " SAME").c_str());
        legend.AddEntry(histograms[i].get(), group.curves[i].legend.c_str(), variable.points ? "p" : "l");
    }
    legend.Draw();
    drawInfo(nEvents);
    save(canvas, outDir + "/" + variable.fileTag + "_" + group.tag);
    return true;
}

bool drawShapes(TFile& file, const Variable& variable, const std::string& yTitle, bool showMean,
                const std::string& nEvents, const std::string& outDir)
{
    struct Species
    {
        std::string positive;
        std::string negative;
        std::string legend;
        Color_t color;
        Style_t lineStyle;
    };
    const std::vector<Species> species = {
        {"pip", "pim", "#pi^{+} + #pi^{-}", kBlack, kSolid},
        {"Kp", "Km", "K^{+} + K^{-}", kGreen + 2, kDashed},
        {"p", "pbar", "p + #bar{p}", kMagenta + 2, kDotted},
    };

    std::vector<std::unique_ptr<TH1D>> histograms;
    std::vector<std::string> legends;
    double maximum = 0.0;
    for (const Species& s : species)
    {
        const TH1D* positive = file.Get<TH1D>((variable.prefix + "_" + s.positive).c_str());
        const TH1D* negative = file.Get<TH1D>((variable.prefix + "_" + s.negative).c_str());
        if (!positive || !negative)
        {
            std::cerr << "No " << variable.prefix << " histogram of " << s.positive << " or " << s.negative
                      << " in " << file.GetName() << "\n";
            return false;
        }
        std::unique_ptr<TH1D> histogram(static_cast<TH1D*>(positive->Clone()));
        histogram->SetDirectory(nullptr);
        histogram->Add(negative);
        legends.push_back(showMean ? s.legend + Form(", #LTp_{T}#GT = %.2f GeV", histogram->GetMean())
                                   : s.legend);
        // divide by the number of particles and by the bin width, so the area is 1
        histogram->Scale(1.0 / histogram->Integral(), "width");
        histogram->SetStats(false);
        histogram->SetLineColor(s.color);
        histogram->SetLineStyle(s.lineStyle);
        histogram->SetLineWidth(2);
        maximum = std::max(maximum, highestBin(*histogram));
        histograms.push_back(std::move(histogram));
    }

    TH1D& frame = *histograms.front();
    frame.SetTitle((variable.title + ", shapes (unit area);" + variable.xTitle + ";" + yTitle).c_str());
    frame.GetYaxis()->SetTitleOffset(1.4);
    frame.SetMinimum(variable.logY ? maximum * 1e-6 : 0.0);
    frame.SetMaximum(variable.logY ? maximum * 50.0 : maximum * 1.5);

    TLegend legend(0.52, 0.68, 0.88, 0.86);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(textSize);
    TCanvas canvas(("c_" + variable.fileTag + "_shapes").c_str(), "", 800, 600);
    prepareCanvas(canvas, variable);
    for (std::size_t i = 0; i < histograms.size(); ++i)
    {
        histograms[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
        legend.AddEntry(histograms[i].get(), legends[i].c_str(), "l");
    }
    legend.Draw();
    drawInfo(nEvents);
    save(canvas, outDir + "/" + variable.fileTag + "_shapes");
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

    const std::string singleDir = std::string(outDir) + "/histograms_single";
    const std::string overlayDir = std::string(outDir) + "/histograms_overlay";
    for (const std::string& dir : {singleDir, overlayDir})
    {
        std::error_code error;
        std::filesystem::create_directories(dir, error);
        if (error)
        {
            std::cerr << "Could not create " << dir << ": " << error.message() << "\n";
            return 1;
        }
    }

    // statistics box: entries, mean, standard deviation, underflow and overflow
    gStyle->SetOptStat("emrou");
    gStyle->SetStatX(0.88);
    gStyle->SetStatY(0.88);
    gStyle->SetStatFormat("6.4g");
    gStyle->SetTitleFontSize(0.045);
    gStyle->SetErrorX(0.0);
    TGaxis::SetMaxDigits(3);

    const TH1D* nCharged = file->Get<TH1D>("hN_ch");
    if (!nCharged)
    {
        std::cerr << "No histogram hN_ch in " << inFileName << "\n";
        return 1;
    }
    const std::string nEvents = std::to_string(std::llround(nCharged->GetEntries()));

    const std::vector<Variable> variables = {
        {"hN", "multiplicity", "Multiplicity", "Multiplicity N", "Events", true, true, false},
        {"hPt", "pt", "Transverse momentum", "p_{T} [GeV]", "Particles / 0.1 GeV", true, false, true},
        {"hEta", "eta", "Pseudorapidity", "#eta", "Particles / 0.1", false, false, true},
        {"hPhi", "phi", "Azimuthal angle", "#phi [rad]", "Particles / 0.098 rad", false, false, true},
    };

    const Color_t positive = kRed + 1;
    const Color_t negative = kAzure + 2;
    const std::vector<Group> groups = {
        {"all and charged particles", "all_ch",
         {
             {"all", "all final-state", "all final-state particles", kBlack, kSolid, kFullCircle},
             {"ch", "charged", "charged particles", kMagenta + 2, kSolid, kFullCircle},
         },
         "charged / all", 0.0, 1.0, true},
        {"#pi^{+} and #pi^{-}", "pions",
         {
             {"pip", "#pi^{+}", "#pi^{+}", positive, kSolid, kFullCircle},
             {"pim", "#pi^{-}", "#pi^{-}", negative, kDashed, kOpenCircle},
         },
         "#pi^{-} / #pi^{+}", 0.9, 1.1, false},
        {"K^{+} and K^{-}", "kaons",
         {
             {"Kp", "K^{+}", "K^{+}", positive, kSolid, kFullCircle},
             {"Km", "K^{-}", "K^{-}", negative, kDashed, kOpenCircle},
         },
         "K^{-} / K^{+}", 0.8, 1.2, false},
        {"p and #bar{p}", "protons",
         {
             {"p", "p", "p", positive, kSolid, kFullCircle},
             {"pbar", "#bar{p}", "#bar{p}", negative, kDashed, kOpenCircle},
         },
         "#bar{p} / p", 0.4, 1.2, false},
    };

    for (const Variable& variable : variables)
    {
        for (const Group& group : groups)
        {
            if (!drawOverlay(*file, variable, group, nEvents, overlayDir))
            {
                return 1;
            }
            for (const Curve& curve : group.curves)
            {
                if (!drawSingle(*file, variable, curve, nEvents, singleDir))
                {
                    return 1;
                }
            }
        }
    }

    const Variable& pt = variables[1];
    const Variable& eta = variables[2];
    if (!drawShapes(*file, pt, "(1/N) dN/dp_{T} [GeV^{-1}]", true, nEvents, overlayDir) ||
        !drawShapes(*file, eta, "(1/N) dN/d#eta", false, nEvents, overlayDir))
    {
        return 1;
    }
    return 0;
}
