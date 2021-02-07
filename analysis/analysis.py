#pylint: disable=too-many-locals,too-many-statements, missing-docstring, pointless-string-statement
from array import array
import yaml
# pylint: disable=import-error, no-name-in-module, unused-import
from ROOT import TH1F, TH2F, TCanvas, TGraph, TLatex, gPad, TFile, TF1
from ROOT import gStyle, gROOT, TStyle, TLegendEntry, TLegend

"""
read predictions from arXiv.1907.12786
"""
def scale(hadron="Omega_ccc", model="SHMC_2021", collision="PbPb", \
        brmode="central"):
    with open(r'prediction.yaml') as fileparam:
        param = yaml.load(fileparam, Loader=yaml.FullLoader)
    yieldmid = param["models"][model][collision][hadron]
    sigma_aa_b = param["statistics"][collision]["sigmaAA_b"]
    lumiaa_monthi_invnb = param["statistics"][collision]["lumiAA_monthi_invnb"]
    nevt = sigma_aa_b * lumiaa_monthi_invnb * 1e9
    bratio = param["branchingratio"][hadron][brmode]
    enhanc = 1
    eff = 1.
    legendtext = '%s N_{ev}(%s) = %.1f B, BR=%.5f%%, #varepsilon=%.2f' \
            % (model, collision, nevt/1e9, bratio*100, eff)
    scale_factor = bratio * enhanc * eff * nevt * yieldmid
    return scale_factor, legendtext

def analysis(hadron="Omega_ccc"):
    gStyle.SetOptStat(0)
    with open(r'prediction.yaml') as fileparam:
        param = yaml.load(fileparam, Loader=yaml.FullLoader)
    models = param["comparison_models"][hadron]["models"]
    collisions = param["comparison_models"][hadron]["collisions"]
    brmode = param["comparison_models"][hadron]["brmode"]
    colors = param["comparison_models"][hadron]["colors"]
    useshape = param["comparison_models"][hadron]["useshape"]
    ymin = param["comparison_models"][hadron]["ymin"]
    ymax = param["comparison_models"][hadron]["ymax"]
    binanal = array('d', param["pt_binning"]["hadron"][hadron])
    dorebin = param["do_corr"]["dorebin"]

    fin = TFile("../Inputs/" + useshape +".root")
    histo_norm = fin.Get("hpred_norm")

    canvas = TCanvas("canvas", "A Simple Graph Example", 881, 176, 668, 616)
    gStyle.SetOptStat(0)
    canvas.SetHighLightColor(2)
    canvas.Range(-1.25, -4.625, 11.25, 11.625)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetBorderSize(2)
    canvas.SetLogy()
    canvas.SetFrameBorderMode(0)
    canvas.SetFrameBorderMode(0)
    canvas.cd()
    gPad.SetLogy()

    hempty = TH2F("hempty", ";p_{T};Yields", 100, 0., 10., 100, ymin, ymax)
    hempty.GetXaxis().SetTitle("p_{T}")
    hempty.GetXaxis().SetLabelFont(42)
    hempty.GetXaxis().SetTitleOffset(1)
    hempty.GetXaxis().SetTitleFont(42)
    hempty.GetYaxis().SetLabelFont(42)
    hempty.GetYaxis().SetTitleOffset(1.35)
    hempty.GetYaxis().SetTitleFont(42)
    hempty.GetZaxis().SetLabelFont(42)
    hempty.GetZaxis().SetTitleOffset(1)
    hempty.GetZaxis().SetTitleFont(42)
    hempty.Draw()
    histolist = [None]*len(models)

    leg = TLegend(0.1471471, 0.6108291, 0.3018018, 0.8747885, "", "brNDC")
    leg.SetBorderSize(1)
    leg.SetLineColor(0)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetTextSize(0.022)
    leg.SetFillStyle(1001)

    for icase, _ in enumerate(models):
        histolist[icase] = histo_norm.Clone("histo_pred%d" % icase)
        scalef, text = scale(hadron, models[icase], collisions[icase], brmode[icase])
        for ibin in range(histolist[icase].GetNbinsX()-1):
            binwdith = histolist[icase].GetBinWidth(ibin+1)
            yvalue = histolist[icase].GetBinContent(ibin+1)
            histolist[icase].SetBinContent(ibin+1, binwdith*scalef*yvalue)
        if dorebin is True:
            histolist[icase] = histolist[icase].Rebin(len(binanal)-1, \
                                "histo_pred%d" % icase, binanal)
        histolist[icase].SetLineColor(colors[icase])
        histolist[icase].SetMarkerColor(colors[icase])
        histolist[icase].SetLineWidth(2)
        histolist[icase].Draw("same")
        text = text + " Yield(tot)=%.2f" % histolist[icase].Integral()
        leg.AddEntry(histolist[icase], text, "pF")
    leg.Draw()
    canvas.SaveAs(hadron+"_results.pdf")
    canvas.SaveAs(hadron+"_results.C")

analysis("Omega_ccc")
analysis("Omega_cc")
analysis("Xi_cc")
analysis("X3872")
