#!/bin/python
########################################################################
# The JTPPPlotterCompPMvsDS by Edgar Kellermann (2016/03/14)
# especially optimised to make comparison plots between trigger and offline plots
# that were created by JTPP
########################################################################

from ROOT import *
import math
import sys
import re
from array import array
import AtlasStyle
#from JTPPPlotter_RatioUtils import divideGraphs

#import JESResponseFitter # Thank you, Antonio!
#import prettyPalette


#ROOT.gROOT.SetBatch(True)

#the 2D plots
  #KEY: TH2D     mjj_vs_BCHCORRCELL_NominalJES;1 Nominal JES
  #KEY: TH2D     pt_vs_BCHCORRCELL_NominalJES;1  Nominal JES

outputfile = TFile("makeResultPdf/JetTriggerPerformancePlots.root","recreate")

## input files
fin1 = TFile.Open("hist-Data16_run300600/hist-Data16_run300600.root")

## sample name in the legend
sampleName1 = "Data 2016, Run 300600"

AtlasStyle.SetAtlasStyle()
gStyle.SetOptStat(0) # get rid of statistics box

# set close or wider Xaxis range
closer = True

# open logfiles from InputHandler to obtain a list of turnons
turnOnData = open("logfile_turnOns.txt") # newest list of JTPP
turnOnList = turnOnData.readlines()

# remove the \n
turnOnList = ([s.strip('\n') for s in turnOnList])
print "Reading logfile_trigger.txt... selected TurnOns are: \n" + str(turnOnList)

# color array if needed
color = [ kRed+2, kBlue+2, kGreen+2, kCyan+2, kMagenta+2, kYellow+2, kRed-2, kBlue+0, kGreen-6, kCyan-2, kMagenta-7, kYellow-7 ]

# define multiHistoList
multiHistoList = ["HLT-j0-1i2c200m8000TLA", "HLT-j0-0i1c500m900TLA"]

# -------------------------- TurnOns ----------------------------------

outputfile.cd()
c = TCanvas("c", sampleName1, 1)

# loop over all multi histograms
for multiHistoName in multiHistoList:

  ## Legend
  legend = TLegend(0.5,0.40,0.85,0.2)
  legend.SetFillColor(0)
  legend.SetTextFont(1)
  legend.SetHeader(sampleName1)
  legend.SetFillStyle(0)


  # multigraph for errors
  MultiGr = TMultiGraph()


  histoCounter = 0

  # Loop over all turnons
  for turnOn in turnOnList:


    ## Decide which turnon will be in which multiHisto
    print "Filling... " + multiHistoName
    if (multiHistoName == "HLT-j0-1i2c200m8000TLA"):
      if (turnOn == "HLT_j0_1i2c200m8000TLA-L1_J100"): print "using... " + turnOn
#      elif (turnOn == "HLT_4j85-HLT_j85"): print "using... " + turnOn
      else: continue

    if (multiHistoName == "HLT-j0-0i1c500m900TLA"):
      if (turnOn == "HLT_j0_0i1c500m900TLA-L1_J100"): print "using... " + turnOn
#      elif (turnOn == "HLT_j110-HLT_j85"): print "using... " + turnOn
      else: continue

    # Loopt over EvaluationType
    for EvalType in ["TDT"]: # ["Emu", "TDT"]:

      histoName1 = "TurnOns" + "/" + "effic_pt_" + EvalType + "_" + turnOn
      histoName1Denum = "TurnOns" + "/" + "effic_DENUM_pt_" + EvalType + "_" + turnOn

      histo1 = fin1.Get(histoName1)
      histo1Denum = fin1.Get(histoName1Denum)

      # Rebin histo
      histo1.Rebin(12)
      histo1Denum.Rebin(12)


      ErrGraph1 = TGraphAsymmErrors(histo1, histo1Denum, 'b(1,1) mode')

      # write in outputfile
      ErrGraph1.Write(turnOn)

      # Divide histograms
      histo1.Divide(histo1, histo1Denum, 1.0, 1.0 , "")

      # get rid of non fully efficient region
      if (turnOn == "HLT_j110-HLT_j85"):
        for n in range (0, 18):
          histo1.SetBinContent(n, -1)
          ErrGraph1.SetPoint(n, 1000, 1000)

      # Set Axis Titles
      if (multiHistoName == "HLT-j0-1i2c200m8000TLA"):
        histo1.GetXaxis().SetTitle("offline m23 [GeV]")
      if (multiHistoName == "HLT-j0-0i1c500m900TLA"):
        histo1.GetXaxis().SetTitle("offline mjj [GeV]")

      histo1.GetYaxis().SetTitle("efficiency")

      # Set Title offset
      histo1.GetXaxis().SetTitleOffset(1.4)
      histo1.GetYaxis().SetTitleOffset(1.4)

      # Set Line style
      histo1.SetLineStyle(0)

      # Set Marker Style
      histo1.SetMarkerStyle(20)
      ErrGraph1.SetMarkerStyle(20)

      # Set Line color
      histo1.SetLineColor(color[histoCounter])
      ErrGraph1.SetLineColor(color[histoCounter])

      # Set marker color
      histo1.SetMarkerColor(color[histoCounter])
      ErrGraph1.SetMarkerColor(color[histoCounter])

      # Set range
      histo1.GetYaxis().SetRangeUser(0.0, 1.1)
      histo1.GetXaxis().SetRangeUser(20.0, 600.0)
      if (closer):
        if (multiHistoName == "HLT-j0-1i2c200m8000TLA"): histo1.GetXaxis().SetRangeUser(20.0, 1000.0)
        if (multiHistoName == "HLT-j0-0i1c500m900TLA"): histo1.GetXaxis().SetRangeUser(20.0, 1200.0)

      # Add legend entry
      legend.AddEntry(histo1,turnOn.rsplit('-',1)[0],"p")

      # check the histoCounter and therefore if Draw need to be with the option "Same" or without
      if (histoCounter == 0):
        histo1.Draw("hist p")
        MultiGr.Add(ErrGraph1)
      else:
        histo1.Draw("Same hist p")
        MultiGr.Add(ErrGraph1)

      histoCounter +=1

      c.Update()

      # Draw legend and ATLAS style label
      legend.Draw("Same")
      AtlasStyle.ATLAS_LABEL(0.5, 0.55, internal = True, preliminary = False, color=1)
      AtlasStyle.myText(0.5, 0.5, 1, "#sqrt{s} = 13 TeV")


      c.Update()

    # Set Title
#    l = TLatex(20,1.11, "#font[12]{" + sampleName1 + "}")
#    l.Draw("Same")

    # Draw multiGraph
    MultiGr.Draw("P")

    c.Update()

# save pdf and eps files
## convert names latex friendly
#triggerPDFNamePrae = re.sub("_","-", turnOn)
#triggerPDFName =  re.sub(".31ETA","-31ETA", triggerPDFNamePrae)

  c.SaveAs("makeResultPdf/plots/"+"efficienciesPT-" + multiHistoName + ".pdf")
  c.SaveAs("makeResultPdf/plots/"+"efficienciesPT-" + multiHistoName + ".pdf")


outputfile.Close()



