#pylint: disable=too-many-locals,too-many-statements, missing-docstring, pointless-string-statement
from math import sqrt
from array import array
import yaml
# pylint: disable=import-error, no-name-in-module, unused-import, too-many-arguments
from ROOT import TH1F, TH2F, TCanvas, TGraph, TLatex, gPad, TFile, TF1
from ROOT import gStyle, gROOT, TStyle, TLegendEntry, TLegend

"""
Macro to perform significance estimation of a given decay channel
in a given collision system. The ingredients are:
    - efficiency vs pt
    - estimated background yield within 3sigma below the peak / event
      bkgperevent is already multiplied by binwidth
    - expected signal yield / event = dN/dpT * binwidth * BR * efficiency
    - to get dN/dpT, I have: Lc cross section from Pythia -> converted to yield by scaling with sigma=70000mub
                              Raa from TAMU and Ncoll
                              -> together they give dN/dpT^AA = dN/dpT^pp * Ncoll * Raa
"""

def analysis(hadron="Lambda_c", collision="pp14p0", yrange="absy1p44", \
        brmode="central", model="Pythia8mode2", centrality="cent010", use_unnorm=1):
    gStyle.SetOptStat(0)
    with open(r'databases/significance.yaml') as filesignificance:
        paramsignificance = yaml.load(filesignificance, Loader=yaml.FullLoader)
    ymin = paramsignificance[hadron][collision][yrange]["ymin"]
    ymax = paramsignificance[hadron][collision][yrange]["ymax"]
    #bin of the final analysis, has to be the binning of efficiency, bkg histos
    binanal = array('d', paramsignificance[hadron][collision][yrange]["binning"])
    nfileyieldth = paramsignificance[hadron][collision][yrange]["theoryfile"]
    nfileraath = paramsignificance[hadron][collision][yrange][centrality]["raafile"]

    nfileeff = paramsignificance[hadron][collision][yrange]["efffile"]
    nhistoeff = paramsignificance[hadron][collision][yrange]["histoeff"]
    nfilebkg = paramsignificance[hadron][collision][yrange][centrality]["bkgfile"]
    nhistobkg = paramsignificance[hadron][collision][yrange][centrality]["histobkg"]
    nhistoyieldth = paramsignificance[hadron][collision][yrange]["histoyield"]
    nhistoyieldth_norm = paramsignificance[hadron][collision][yrange]["histoyield_norm"]
    nhistoraath = paramsignificance[hadron][collision][yrange][centrality]["historaa"]

    with open(r'databases/general.yaml') as fileparamgen:
        paramgen = yaml.load(fileparamgen, Loader=yaml.FullLoader)
    with open(r'databases/theory_yields.yaml') as fileyields:
        paramyields = yaml.load(fileyields, Loader=yaml.FullLoader)

    textcollision = paramgen["text_string"][collision]
    textrapid = paramgen["text_string"][yrange]
    textmodel = paramgen["text_string"][model]

    sigma_aa_b = paramgen["statistics"][collision]["sigmaAA_b"]
    lumiaa_monthi_invnb = paramgen["statistics"][collision]["lumiAA_monthi_invnb"]
    event_scaling = paramgen["statistics"][collision][centrality]["EventScaling"]
    #nevt = sigma_aa_b * 5.6 * 1e9 / 10 # Sourav: 5.6 nb for one year of run 5, assume constant centrality dist. -> scale with 10 to have 0-10%
    nevt = sigma_aa_b * lumiaa_monthi_invnb * 1e9 * 7 / event_scaling # 5*7=35 nb for whole Run 5, assume constant centrality dist. -> scale with 5 to have 30-50%
    #nevt = 2.*1e9
    bratio = paramgen["branchingratio"][hadron][brmode]
    decaychannel = paramgen["latexparticle"][hadron]
    ncoll = paramgen["statistics"][collision][centrality]["Ncoll"]

    #yieldmid = paramyields[model][collision][yrange][hadron]
    #note: here I hardcode "absy0p5" because in the database file we only hvae TAMU theory calculations for this rapidity
    yieldmid = paramyields[model][collision]["absy0p5"][hadron]
    y_corr = paramyields[model][collision]["absy0p5"]["y_corr"] # correction for TAMU which is in |y|<0.5
    y_corrPythia = paramsignificance[hadron][collision][yrange]["y_corrPythia"] #correction factor for pp yield which is |y|<1.44  
    #y_corrPythia = 1.0
    #y_corrPythia = 1.0/1.44 #correction factor for pp yield which is |y|<1.44
    energy_corr = 5.02/14 #correction factor because pp yield is @14 TeV and PbPb is @5 TeV
    text = '%s, N_{ev} = %.0f 10^{12}' % (textmodel, nevt/1e12)
    text_a = '%s, %s, BR=%.2f%%' % (decaychannel, textrapid, bratio*100)
    text_b = 'ALICE3 projection, with IRIS, no PID, %s' % textcollision
    fileeff = TFile(nfileeff)
    histoeff = fileeff.Get(nhistoeff)
    filebkg = TFile(nfilebkg)
    hbkgperevent = filebkg.Get(nhistobkg)

    fileyieldth = TFile(nfileyieldth)
    histoyieldth = None

    if use_unnorm == 1:
        histodndptth = fileyieldth.Get(nhistoyieldth)
        histodndptth.Scale(1./70000.) #TEMPORARY this is a fix to account for the
                                      #conversion from a cross-section in mub
                                      #to yields, sigma=70000 mub
    else:
        histodndptth = fileyieldth.Get(nhistoyieldth_norm)
        histodndptth.Scale(yieldmid*y_corr)

    histoyieldth = histodndptth.Clone("histoyieldth")

    fileraath = TFile(nfileraath);
    historaath = fileraath.Get(nhistoraath)

    deltay = paramsignificance[hadron][collision][yrange]["deltay"] 
    for ibin in range(histoyieldth.GetNbinsX()):
        binwdith = histoyieldth.GetBinWidth(ibin+1)*deltay
        yieldperevent = histoyieldth.GetBinContent(ibin+1)*y_corrPythia*energy_corr*binwdith*bratio*y_corr*ncoll*historaath.GetBinContent(ibin+1)
        histoyieldth.SetBinContent(ibin+1, yieldperevent)
        histoyieldth.SetBinError(ibin+1, 0.)
    histoyieldth = histoyieldth.Rebin(len(binanal)-1, \
             "histoyieldth", binanal)
    histosignfperevent = histoyieldth.Clone("histosignfperevent")
    histosignf = histoyieldth.Clone("histosignf")
    histosignal = histoyieldth.Clone("histosignal")
    histosigoverbkg = histoyieldth.Clone("histosigoverbkg")

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

    hempty = TH2F("hempty", ";p_{T} (GeV/c); Significance(3#sigma)", 100, 0., 10., 100, ymin, ymax)
    hempty.GetXaxis().SetTitle("p_{T} (GeV/c)")
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

    histosignf = histosignfperevent.Clone("histosignf")
    for ibin in range(histoyieldth.GetNbinsX()):
        yieldperevent = histoyieldth.GetBinContent(ibin+1)
        bkgperevent = hbkgperevent.GetBinContent(ibin+1)
        eff = histoeff.GetBinContent(ibin+1)
        #signalperevent = yieldperevent
        signalperevent = eff*yieldperevent
        significanceperevent = 0
        if bkgperevent > 0:
            significanceperevent = signalperevent/sqrt(signalperevent+bkgperevent)
        signaloverbkg = 0
        if bkgperevent > 0:
            signaloverbkg = signalperevent/bkgperevent
        histosignfperevent.SetBinContent(ibin+1, significanceperevent)
        histosignfperevent.SetBinError(ibin+1, 0.)
        histosignf.SetBinContent(ibin+1, significanceperevent*sqrt(nevt))
        histosignf.SetBinError(ibin+1, 0.)
        histosignal.SetBinContent(ibin+1, signalperevent*nevt)
        if histosignf.GetBinContent(ibin+1) != 0:
            histosignal.SetBinError(ibin+1, 1./histosignfperevent.GetBinContent(ibin+1))
        else:
            histosignal.SetBinError(ibin+1, 0.)
        histosigoverbkg.SetBinContent(ibin+1, signaloverbkg)
        histosigoverbkg.SetBinError(ibin+1, 0.)


    histosignfperevent.SetLineColor(1)
    histosignfperevent.SetMarkerColor(1)
    histosignfperevent.SetLineWidth(1)
    histosignf.SetLineColor(1)
    histosignf.SetMarkerColor(1)
    histosignf.SetLineWidth(2)
    histosignf.Draw("same")
    t_b = TLatex()
    t_b.SetNDC()
    t_b.SetTextFont(42)
    t_b.SetTextColor(1)
    t_b.SetTextSize(0.035)
    t_b.SetTextAlign(12)
    t_b.DrawLatex(0.2, 0.85, text_b)
    t_c = TLatex()
    t_c.SetNDC()
    t_c.SetTextFont(42)
    t_c.SetTextColor(1)
    t_c.SetTextSize(0.03)
    t_c.SetTextAlign(12)
    t_c.DrawLatex(0.2, 0.80, text)
    t_a = TLatex()
    t_a.SetNDC()
    t_a.SetTextFont(42)
    t_a.SetTextColor(1)
    t_a.SetTextSize(0.035)
    t_a.SetTextAlign(12)
    t_a.DrawLatex(0.2, 0.75, text_a)
    canvas.SaveAs(hadron+"_results.pdf")
    canvas.SaveAs(hadron+"_results.C")

    foutput = TFile("scopingdoc/foutput" + hadron + centrality + "_d9n_y1p44.root", "recreate")
    #foutput = TFile("foutput" + hadron + centrality + "_d9n_y1p0_test.root", "recreate")
    foutput.cd()
    histoeff.Write()
    hbkgperevent.Write()
    histosignfperevent.Write()
    histoyieldth.Write()
    histosignf.Write()
    histosignal.Write()
    histodndptth.Write()
    histosigoverbkg.Write()
#analysis("Lambda_c", "PbPb5p02", "absy1p44", "central", "TAMU", "cent010", 1)
analysis("Lambda_c", "PbPb5p02", "absy1p44", "central", "TAMU", "cent3050", 1)

#analysis("Lambda_c", "PbPb5p02", "absy1p0", "central", "TAMU", "cent010", 1)
#analysis("Lambda_c", "PbPb5p02", "absy1p0", "central", "TAMU", "cent3050", 1)
#analysis("Lambda_c", "PbPb2p76", "absy1p0", "central", "Stat_ChoLee_2", 0)
