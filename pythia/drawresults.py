#!/usr/bin/env python
import os

from ROOT import TF1, TH2F, TCanvas, TFile, TLatex, gROOT, gStyle


def makeSavePaths(title, *fileFormats, outputdir="outputPlots"):
    """
    Saves the canvas as the desired output format in an output directory (default = outputPlots)
    """
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    return [outputdir + "/" + title + fileFormat for fileFormat in fileFormats]


def plotsinglevar(
    filename="tot.root",
    histonmasssig="hInvmassvsGluonE",
    xminproj=20.,
    xmaxproj=40.,
    xmin=0.,
    xmax=30.,
    rebinx=5,
    rebiny=1,
    logx=0,
    logy=0,
    xaxis="m_{c#bar{c}} (GeV)",
    latex="E_{g}",
):
    gStyle.SetOptStat(0)
    gROOT.SetBatch(1)

    fileSig = TFile(filename)
    histo2d = fileSig.Get("%s" % (histonmasssig))
    print("NbinsX=", histo2d.GetXaxis().GetNbins())
    print("NbinsY=", histo2d.GetYaxis().GetNbins())
    minbin = histo2d.GetXaxis().FindBin(xminproj)
    maxbin = histo2d.GetXaxis().FindBin(xmaxproj)
    hvar = histo2d.ProjectionY("hvar", minbin, maxbin)
    hvar.GetXaxis().SetRangeUser(xmin, xmax)
    hvar.Rebin(rebinx);
    hvar.Draw()
    ymax = hvar.GetMaximum()*1.5;
    hempty = TH2F("hempty", ";%s; Entries" % xaxis, 100, xmin, xmax, 100, 0.1, ymax)
    hempty.GetXaxis().SetLabelFont(42)
    hempty.GetXaxis().SetTitleOffset(1)
    hempty.GetXaxis().SetLabelSize(0.03)
    hempty.GetXaxis().SetTitleFont(42)
    hempty.GetYaxis().SetLabelFont(42)
    hempty.GetYaxis().SetTitleOffset(1.35)
    hempty.GetYaxis().SetTitleFont(42)
    hempty.GetZaxis().SetLabelFont(42)
    hempty.GetZaxis().SetTitleOffset(1)
    hempty.GetZaxis().SetTitleFont(42)

    canvas = TCanvas("canvas", "A Simple Graph Example", 881, 176, 668, 616)
    gStyle.SetOptStat(0)
    canvas.SetHighLightColor(2)
    canvas.Range(-1.25, -4.625, 11.25, 11.625)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetBorderSize(2)
    canvas.SetFrameBorderMode(0)
    canvas.SetFrameBorderMode(0)
    if logx == 1:
        canvas.SetLogx()
    if logy == 1:
        canvas.SetLogy()
    canvas.cd()
    hempty.Draw("")
    hvar.Draw("PEsame")
    latexa = TLatex()
    latexa.SetTextSize(0.04)
    latexa.SetTextFont(42)
    latexa.SetTextAlign(3)
    xave = xmin + (xmax - xmin) / 4.0
    latexa.DrawLatex(
        xave, ymax * 0.7, "%.1f < p_{T} (%s) < %.1f GeV" % (xminproj, latex, xmaxproj)
    )
    canvas.SaveAs("%s%.2f_%.2f.pdf" % (histonmasssig, xminproj, xmaxproj))

plotsinglevar("tot.root", "hInvmassvsGluonE", 5., 10., 0., 30., 5, 1, 0, 0,
              "m_{c#bar{c}} (GeV)", "E_{g}")
plotsinglevar("tot.root", "hInvmassvsGluonE", 10., 20., 0., 30., 5, 1, 0, 0,
              "m_{c#bar{c}} (GeV)", "E_{g}")
plotsinglevar("tot.root", "hInvmassvsGluonE", 20., 40., 0., 30., 5, 1, 0, 0,
              "m_{c#bar{c}} (GeV)", "E_{g}")
plotsinglevar("tot.root", "hInvmassvsGluonE", 0., 9999., 0., 30., 5, 1, 0, 0,
              "m_{c#bar{c}} (GeV)", "E_{g}")
