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
fin1 = TFile.Open("hist-NewLUTTestsROI7601/hist-NewLUTTestsROI7601.root") 
fin2 = TFile.Open("hist-NewLUTTestsROI7602/hist-NewLUTTestsROI7602.root") 

## sample name in the legend and title
title = "LUT"
sampleName1 = "non linear"
sampleName2 = "linear"

AtlasStyle.SetAtlasStyle()
gStyle.SetOptStat(0) # get rid of statistics box

# set close or wider Xaxis range
closer = True

# open logfiles from InputHandler to see what we need to read out
turnOnData = open("logfile_turnOns_LUT.txt") # newest list of JTPP
#turnOnData = open("logfile_turnOns.txt") # newest list of JTPP
turnOnList = turnOnData.readlines()

# remove the \n
turnOnList = ([s.strip('\n') for s in turnOnList])
print "Reading logfile_trigger.txt... selected TurnOns are: \n" + str(turnOnList)

# color array if needed
color = [ kRed+2, kBlue+2, kGreen+2, kCyan+2, kMagenta+2, kYellow+2, kRed-2, kBlue+0, kGreen-6, kCyan-2, kMagenta-7, kYellow-7 ]


# -------------------------- TurnOns ----------------------------------

for turnOn in turnOnList:  # triggerList: WE ONLY LOOK AT J60

  ### mapping procedure

  turnOnMapped = turnOn
  if (turnOn == "L1_3J40-HLT_j60"): turnOnMapped = "L1_3J70-HLT_j60"
  if (turnOn == "L1_4J15-HLT_j45"): turnOnMapped = "L1_4J37-HLT_j45"
  if (turnOn == "L1_J100-HLT_j85"): turnOnMapped = "L1_J146-HLT_j85"
  if (turnOn == "L1_J15-HLT_j15"): turnOnMapped = "L1_J37-HLT_j15"
  if (turnOn == "L1_J20-HLT_j25"): turnOnMapped = "L1_J44-HLT_j25"
  if (turnOn == "L1_J75-HLT_j45"): turnOnMapped = "L1_J114-HLT_j45"
  if (turnOn == "L1_J75-HLT_j60"): turnOnMapped = "L1_J114-HLT_j60"
  if (turnOn == "L1_J15.31ETA49-HLT_j15_320eta490"): turnOnMapped = "L1_J37.31ETA49-HLT_j15_320eta490"
  if (turnOn == "L1_J40.31ETA49-HLT_j15_320eta490"): turnOnMapped = "L1_J70.31ETA49-HLT_j15_320eta490"

  print "turnOn:" + turnOn
  print "turnOnMapped:" + turnOnMapped

  outputfile.cd()
  c = TCanvas()

  ## Legend
  legend = TLegend(0.5,0.40,0.75,0.25)
  legend.SetFillColor(0)
  legend.SetTextFont(1)
  legendName1 = turnOn.split("-")
  legendName2 = turnOnMapped.split("-")
  legend.SetHeader(legendName1[0] + " (" + legendName2[0] + ")")
  legend.SetFillStyle(0)

  # multigraph for errors
  MultiGr = TMultiGraph()

  # Loopt over EvaluationType
  for EvalType in ["Emu", "TDT"]:

    histoName1 = "TurnOns" + "/" + "effic_pt_" + EvalType + "_" + turnOnMapped #MODIFIED because of mapping
    histoName1Denum = "TurnOns" + "/" + "effic_DENUM_pt_" + EvalType + "_" + turnOnMapped #MODIFIED because of mapping
    histoName2 = "TurnOns" + "/" + "effic_pt_" + EvalType + "_" + turnOn
    histoName2Denum = "TurnOns" + "/" + "effic_DENUM_pt_" + EvalType + "_" + turnOn

    histo1 = fin1.Get(histoName1)
    histo1Denum = fin1.Get(histoName1Denum)
    histo2 = fin2.Get(histoName2)
    histo2Denum = fin2.Get(histoName2Denum)

    ErrGraph1 = TGraphAsymmErrors(histo1, histo1Denum, 'b(1,1) mode')
    ErrGraph2 = TGraphAsymmErrors(histo2, histo2Denum, 'b(1,1) mode')

    # Divide histograms
    histo1.Divide(histo1, histo1Denum, 1.0, 1.0 , "")
    histo2.Divide(histo2, histo2Denum, 1.0, 1.0 , "")

    # Set Axis Titles
    histo1.GetXaxis().SetTitle("offline pt [GeV]")
    histo1.GetYaxis().SetTitle("efficiency")

    # Set Title offset
    histo1.GetXaxis().SetTitleOffset(1.4)
    histo1.GetYaxis().SetTitleOffset(1.4)

    histo2.GetXaxis().SetTitleOffset(1.4)
    histo2.GetYaxis().SetTitleOffset(1.4)

    # Set Line style
    histo1.SetLineStyle(0)
    histo2.SetLineStyle(0)

    # Set Marker Style
    if (EvalType == "Emu"):
      histo1.SetMarkerStyle(20)
      ErrGraph1.SetMarkerStyle(20)
    else:
      histo1.SetMarkerStyle(24)
      ErrGraph1.SetMarkerStyle(24)

    if (EvalType == "Emu"):
      histo2.SetMarkerStyle(20)
      ErrGraph2.SetMarkerStyle(20)
    else:
      histo2.SetMarkerStyle(24)
      ErrGraph2.SetMarkerStyle(24)

    # Set Line color
    histo1.SetLineColor(kRed+2)
    ErrGraph1.SetLineColor(kRed+2)
    histo2.SetLineColor(kBlue+2)
    ErrGraph2.SetLineColor(kBlue+2)

    # Set marker color
    histo1.SetMarkerColor(kRed+2)
    ErrGraph1.SetMarkerColor(kRed+2)
    histo2.SetMarkerColor(kBlue+2)
    ErrGraph2.SetMarkerColor(kBlue+2)

    # Set range
    histo1.GetYaxis().SetRangeUser(0.0, 1.1)
    histo1.GetXaxis().SetRangeUser(20.0, 300.0)
    if (closer):
      if (turnOn == "HLT_j25-HLT_j15"): histo1.GetXaxis().SetRangeUser(20.0, 200.0)
      if (turnOn == "HLT_j60-HLT_j25"): histo1.GetXaxis().SetRangeUser(20.0, 200.0)
      if (turnOn == "HLT_j110-HLT_j85"): histo1.GetXaxis().SetRangeUser(20.0, 250.0)
      if (turnOn == "HLT_j360-HLT_j260"): histo1.GetXaxis().SetRangeUser(20.0, 600.0)
      if (turnOn == "L1_J15-HLT_j15"): histo1.GetXaxis().SetRangeUser(20.0, 200.0)
      if (turnOn == "L1_J20-HLT_j25"): histo1.GetXaxis().SetRangeUser(20.0, 200.0)
      if (turnOn == "L1_J100-HLT_j85"): histo1.GetXaxis().SetRangeUser(60.0, 300.0)
      if (turnOn == "HLT_3j175-HLT_j110"): histo1.GetXaxis().SetRangeUser(20.0, 350.0)
      if (turnOn == "HLT_4j85-HLT_j85"): histo1.GetXaxis().SetRangeUser(20.0, 250.0)
      if (turnOn == "HLT_6j45-HLT_j60"): histo1.GetXaxis().SetRangeUser(20.0, 200.0)
      if (turnOn == "L1_3J40-HLT_j60"): histo1.GetXaxis().SetRangeUser(20.0, 200.0)
      if (turnOn == "L1_4J15-HLT_j45"): histo1.GetXaxis().SetRangeUser(20.0, 200.0)
      if (turnOn == "HLT_j25_320eta490-HLT_j15_320eta490"): histo1.GetXaxis().SetRangeUser(20.0, 250.0)
      if (turnOn == "HLT_j60_320eta490-HLT_j25_320eta490"): histo1.GetXaxis().SetRangeUser(20.0, 250.0)
      if (turnOn == "L1_J15.31ETA49-HLT_j15_320eta490"): histo1.GetXaxis().SetRangeUser(20.0, 200.0)
      if (turnOn == "L1_J40.31ETA49-HLT_j15_320eta490"): histo1.GetXaxis().SetRangeUser(20.0, 200.0)


    if (EvalType == "Emu"):
      legend.AddEntry(histo1, sampleName1 + " (" + EvalType + ")","lp")


    # Drawing (Histos) and Adding (Graphs)
    if (EvalType == "Emu"):
      histo1.Draw("hist p")
      MultiGr.Add(ErrGraph1)
#    else:
#      histo1.Draw("Same p")
#      MultiGr.Add(ErrGraph1) 
    if (("31ETA49" in turnOn)&(EvalType == "Emu")):
      histo2.Draw("same hist p") #don't draw TDT here for mapped one
      MultiGr.Add(ErrGraph2)
      legend.AddEntry(histo2,sampleName2 + " (" + EvalType + ")","lp")
    if ("31ETA49" not in turnOn):
      histo2.Draw("same ghist p")
      MultiGr.Add(ErrGraph2)
      legend.AddEntry(histo2,sampleName2 + " (" + EvalType + ")","lp")


    c.Update()

    # Draw legend and ATLAS style label
    legend.Draw("Same")
    AtlasStyle.ATLAS_LABEL(0.5, 0.55, internal = False, preliminary = True, color=1)
    AtlasStyle.myText(0.5, 0.5, 1, "#sqrt{s} = 13 TeV")


    c.Update()

  # Set Title
  l = TLatex(20, 1.11, "#font[12]{" + title + "}")
  l.Draw("Same")

  # Draw multiGraph
  MultiGr.Draw("P")

  c.Update()

  # save pdf and eps files
  ## convert names latex friendly
  triggerPDFNamePrae = re.sub("_","-", turnOn)
  triggerPDFName =  re.sub(".31ETA","-31ETA", triggerPDFNamePrae)

  c.SaveAs("makeResultPdf/plots/"+"efficienciesPT-" + EvalType + "-" + triggerPDFName + ".pdf")
  c.SaveAs("makeResultPdf/plots/"+"efficienciesPT-" + EvalType + "-" + triggerPDFName + ".pdf")


outputfile.Close()



