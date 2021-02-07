#pylint: disable=too-many-locals,too-many-statements, missing-docstring, pointless-string-statement
from array import array
import yaml
import pandas as pd #pylint: disable=import-error
# pylint: disable=import-error, no-name-in-module, unused-import
from ROOT import TH1F, TH2F, TCanvas, TGraph, TLatex, gPad, TFile, TF1
from ROOT import gStyle, gROOT, TStyle, TLegendEntry

"""
read predictions from arXiv.1907.12786
"""

def read_predictions(hadron="Omega_ccc", collision="PbPb"):
    do2pipnorm = 1
    with open(r'prediction.yaml') as fileparam:
        param = yaml.load(fileparam, Loader=yaml.FullLoader)
    latexname = param["latexparticle"][hadron]
    filename = param["ptshape"][hadron]["inputfile"]
    miny = param["ptshape"][hadron]["miny"]
    maxy = param["ptshape"][hadron]["maxy"]
    nbins = param["pt_binning"]["nbins"]
    minb = param["pt_binning"]["minb"]
    maxb = param["pt_binning"]["maxb"]
    width = (maxb-minb)/float(nbins)
    energy = 2.76

    legendtext = '%s %.2f TeV, arXiv.1907.12786 (Coal.2)' % (collision, energy)
    dfpred = pd.read_csv("../Inputs/" + filename+".csv")
    ptcand = dfpred["pt"]
    cross = dfpred["cross"]

    ptcand_val = array('f', ptcand)
    cross_val = array('f', cross)
    if do2pipnorm == 1:
        cross_val_ptscaled = array('f', [2*3.14*a*b for a, b in zip(ptcand_val, cross_val)])
    if do2pipnorm == 0:
        cross_val_ptscaled = cross_val
    npoints = 25
    grpred = TGraph(npoints, ptcand_val, cross_val_ptscaled)
    grpred.SetLineColor(1)
    grpred.SetLineWidth(4)
    grpred.SetMarkerColor(1)
    grpred.SetMarkerStyle(21)
    grpred.SetTitle('')
    grpred.GetXaxis().SetTitle('p_{T} (GeV)')
    grpred.GetYaxis().SetTitle('dN/dp_{T} (%s) (GeV^{-1})' % latexname)

    histo = TH1F("hpred", ";p_{T}; #Delta N/#Delta p_{T}", nbins, minb, maxb)
    norm = 0.
    for i in range(nbins-1):
        histo.SetBinContent(i+1, grpred.Eval(i*width+width/2.))
        print(i+1, i*width+width/2., grpred.Eval(i*width+width/2.))
        norm = norm + width*grpred.Eval(i*width+width/2.)
    fout = TFile("../Inputs/" + hadron+".root", "recreate")
    fout.cd()
    grpred.Write()
    histo.Write()
    histo_norm = TH1F("hpred_norm", ";p_{T}; #Delta N/#Delta p_{T}", nbins, minb, maxb)
    for i in range(nbins-1):
        histo_norm.SetBinContent(i+1, histo.GetBinContent(i+1)/norm)
    histo_norm.Write()

    cpred = TCanvas("cpred", "A Simple Graph Example", 881, 176, 668, 616)
    gStyle.SetOptStat(0)
    cpred.SetHighLightColor(2)
    cpred.Range(-1.25, -4.625, 11.25, 11.625)
    cpred.SetFillColor(0)
    cpred.SetBorderMode(0)
    cpred.SetBorderSize(2)
    cpred.SetLogy()
    cpred.SetFrameBorderMode(0)
    cpred.SetFrameBorderMode(0)
    cpred.cd()
    gPad.SetLogy()
    hempty = TH2F("hempty", ";p_{T};dN/d p_{T}", 100, 0., 10., 100, miny, maxy)
    hempty.GetYaxis().SetTitle('dN/dp_{T} (%s) (GeV^{-1})' % latexname)
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
    grpred.Draw('Psame')
    histo.SetLineColor(1)
    histo.SetLineWidth(2)
    histo.Draw("same")
    tex = TLatex()
    tex.SetNDC()
    tex.SetTextFont(40)
    tex.SetTextColor(1)
    tex.SetTextSize(0.03)
    tex.SetTextAlign(12)
    tex.DrawLatex(0.12, 0.85, legendtext + ", dNdy=%.8f" % norm)
    #tex.Draw()
    cpred.SaveAs("../Inputs/" + filename+".pdf")

read_predictions("Omega_ccc", "PbPb")
read_predictions("Omega_cc", "PbPb")
read_predictions("Xi_cc", "PbPb")
read_predictions("X3872", "PbPb")
read_predictions("Lambda_c", "PbPb")
