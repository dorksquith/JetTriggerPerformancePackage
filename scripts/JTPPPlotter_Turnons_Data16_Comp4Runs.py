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

outputfile = TFile("makeResultPdf/Turnons.root","recreate")

## input files
fin1 = TFile.Open("hist-Data16_run297730/hist-Data16_run297730.root")
fin2 = TFile.Open("hist-Data16_run298967/hist-Data16_run298967.root")
fin3 = TFile.Open("hist-Data16_run300600/hist-Data16_run300600.root")
fin4 = TFile.Open("hist-Data16_run300800/hist-Data16_run300800.root")


sampleName1 = "49"
sampleName2 = "589"
sampleName3 = "1812"
sampleName4 = "2028"


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
#multiHistoList = ["HLT_3j175-HLT_j110", "HLT_4j85-HLT_j85", "HLT_6j45-HLT_j60", "HLT_j25-HLT_j15", "HLT_j60-HLT_j25", "HLT_j110-HLT_j85", "L1_J15-HLT_j15", "L1_J50-HLT_j25", "L1_J75-HLT_j45", "L1_J100-HLT_j85", "L1_J20.31ETA49-HLT_j15_320eta490", "L1_J50.31ETA49-HLT_j45_320eta490", "L1_J75.31ETA49-HLT_j60_320eta490", "L1_3J40-HLT_j60", "L1_4J15-HLT_j45"]

multiHistoList = ["HLT_j25-HLT_j15", "HLT_j60-HLT_j25", "HLT_j110-HLT_j85", "L1_J15-HLT_j15", "L1_J50-HLT_j25", "L1_J75-HLT_j45", "L1_J100-HLT_j85", "L1_J75.31ETA49-HLT_j60_320eta490", "L1_J75.31ETA49-HLT_j45_320eta490", "L1_J20.31ETA49-HLT_j15_320eta490", "L1_J50.31ETA49-HLT_j45_320eta490"]

# -------------------------- TurnOns ----------------------------------

outputfile.cd()
c = TCanvas("c", sampleName1, 1)

# loop over all multi histograms
for multiHistoName in multiHistoList:

  ## Legend
  legend = TLegend(0.45,0.50,0.9,0.2)
  legend.SetFillColor(0)
  legend.SetTextFont(1)
  legend.SetHeader("# of bunches")
  legend.SetFillStyle(0)


  # multigraph for errors
  MultiGr = TMultiGraph()


  histoCounter = 0

  # Loop over all turnons
  for turnOn in turnOnList:

    print "Filling... " + multiHistoName

    if (multiHistoName == turnOn) : print "using... " + turnOn
    else: continue

#    if (multiHistoName == "HLT_3j175-HLT_j110"):
#      if (turnOn == "HLT_3j175-HLT_j110"): print "using... " + turnOn
#      elif (turnOn == "HLT_4j85-HLT_j85"): print "using... " + turnOn
#      elif (turnOn == "HLT_6j45-HLT_j60"): print "using... " + turnOn
#      else: continue

#    if (multiHistoName == "HLT-singleJetTriggers"):
#      if (turnOn == "HLT_j25-HLT_j15"): print "using... " + turnOn
#      elif (turnOn == "HLT_j60-HLT_j25"): print "using... " + turnOn
#      elif (turnOn == "HLT_j110-HLT_j85"): print "using... " + turnOn
#      else: continue
#
#    if (multiHistoName == "L1-singleJetTriggers"):
#      if (turnOn == "L1_J15-HLT_j15"): print "using... " + turnOn
#      elif (turnOn == "L1_J50-HLT_j25"): print "using... " + turnOn
#      elif (turnOn == "L1_J75-HLT_j45"): print "using... " + turnOn
#      elif (turnOn == "L1_J100-HLT_j85"): print "using... " + turnOn
#      else: continue
#
#    if (multiHistoName == "L1-forwardJetTriggers"):
#      if (turnOn == "L1_J20.31ETA49-HLT_j15_320eta490"): print "using... " + turnOn
#      elif (turnOn == "L1_J50.31ETA49-HLT_j45_320eta490"): print "using... " + turnOn
#      elif (turnOn == "L1_J75.31ETA49-HLT_j60_320eta490"): print "using... " + turnOn
#      else: continue
#
#    if (multiHistoName == "L1-multiJetTriggers"):
#      if (turnOn == "L1_3J40-HLT_j60"): print "using... " + turnOn
#      elif (turnOn == "L1_4J15-HLT_j45"): print "using... " + turnOn
#      else: continue


    # Loopt over EvaluationType
    for EvalType in ["Emu"]: # ["Emu", "TDT"]:

      histoName1 = "TurnOns" + "/" + "effic_pt_" + EvalType + "_" + turnOn 
      histoName1Denum = "TurnOns" + "/" + "effic_DENUM_pt_" + EvalType + "_" + turnOn

      histo1 = fin1.Get(histoName1)
      histo1Denum = fin1.Get(histoName1Denum)

      histo2 = fin2.Get(histoName1)
      histo2Denum = fin2.Get(histoName1Denum)

      histo3 = fin3.Get(histoName1)
      histo3Denum = fin3.Get(histoName1Denum)

      histo4 = fin4.Get(histoName1)
      histo4Denum = fin4.Get(histoName1Denum)

      # Rebin
      histo1.Rebin(4)
      histo1Denum.Rebin(4)
      histo2.Rebin(4)
      histo2Denum.Rebin(4)
      histo3.Rebin(4)
      histo3Denum.Rebin(4)
      histo4.Rebin(4)
      histo4Denum.Rebin(4)

      ErrGraph1 = TGraphAsymmErrors(histo1, histo1Denum, 'b(1,1) mode')
      ErrGraph2 = TGraphAsymmErrors(histo2, histo2Denum, 'b(1,1) mode')
      ErrGraph3 = TGraphAsymmErrors(histo3, histo3Denum, 'b(1,1) mode')
      ErrGraph4 = TGraphAsymmErrors(histo4, histo4Denum, 'b(1,1) mode')

      # Divide histograms
      histo1.Divide(histo1, histo1Denum, 1.0, 1.0 , "")
      histo2.Divide(histo2, histo2Denum, 1.0, 1.0 , "")
      histo3.Divide(histo3, histo3Denum, 1.0, 1.0 , "")
      histo4.Divide(histo4, histo4Denum, 1.0, 1.0 , "")

      # get rid of non fully efficient region
      if (turnOn == "HLT_j110-HLT_j85"):
        for n in range (0, 16):
          histo1.SetBinContent(n, -1)
          ErrGraph1.SetPoint(n, 1000, 1000)

          histo2.SetBinContent(n, -1)
          ErrGraph2.SetPoint(n, 1000, 1000)

          histo3.SetBinContent(n, -1)
          ErrGraph3.SetPoint(n, 1000, 1000)

          histo4.SetBinContent(n, -1)
          ErrGraph4.SetPoint(n, 1000, 1000)

      if (turnOn == "L1_J75.31ETA49-HLT_j60_320eta490"):
        for n in range (0, 16):
          histo1.SetBinContent(n, -1)
          ErrGraph1.SetPoint(n, 1000, 1000)

          histo2.SetBinContent(n, -1)
          ErrGraph2.SetPoint(n, 1000, 1000)

          histo3.SetBinContent(n, -1)
          ErrGraph3.SetPoint(n, 1000, 1000)

          histo4.SetBinContent(n, -1)
          ErrGraph4.SetPoint(n, 1000, 1000)

      if (turnOn == "L1_J50.31ETA49-HLT_j45_320eta490"):
        for n in range (0, 13):
          histo1.SetBinContent(n, -1)
          ErrGraph1.SetPoint(n, 1000, 1000)

          histo2.SetBinContent(n, -1)
          ErrGraph2.SetPoint(n, 1000, 1000)

          histo3.SetBinContent(n, -1)
          ErrGraph3.SetPoint(n, 1000, 1000)

          histo4.SetBinContent(n, -1)
          ErrGraph4.SetPoint(n, 1000, 1000)

      if (turnOn == "HLT_j60-HLT_j25"):
        for n in range (0, 7):
          histo1.SetBinContent(n, -1)
          ErrGraph1.SetPoint(n, 1000, 1000)

          histo2.SetBinContent(n, -1)
          ErrGraph2.SetPoint(n, 1000, 1000)

          histo3.SetBinContent(n, -1)
          ErrGraph3.SetPoint(n, 1000, 1000)

          histo4.SetBinContent(n, -1)
          ErrGraph4.SetPoint(n, 1000, 1000)

      if (turnOn == "L1_J75-HLT_j45"):
        for n in range (0, 15):
          histo1.SetBinContent(n, -1)
          ErrGraph1.SetPoint(n, 1000, 1000)

          histo2.SetBinContent(n, -1)
          ErrGraph2.SetPoint(n, 1000, 1000)

          histo3.SetBinContent(n, -1)
          ErrGraph3.SetPoint(n, 1000, 1000)

          histo4.SetBinContent(n, -1)
          ErrGraph4.SetPoint(n, 1000, 1000)

      if (turnOn == "L1_J100-HLT_j85"):
        for n in range (0, 15):
          histo1.SetBinContent(n, -1)
          ErrGraph1.SetPoint(n, 1000, 1000)

          histo2.SetBinContent(n, -1)
          ErrGraph2.SetPoint(n, 1000, 1000)

          histo3.SetBinContent(n, -1)
          ErrGraph3.SetPoint(n, 1000, 1000)

          histo4.SetBinContent(n, -1)
          ErrGraph4.SetPoint(n, 1000, 1000)


      # Set Axis Titles
      histo1.GetXaxis().SetTitle("offline pt [GeV]")
      histo1.GetYaxis().SetTitle("efficiency")

      # Set Title offset
      histo1.GetXaxis().SetTitleOffset(1.4)
      histo1.GetYaxis().SetTitleOffset(1.4)

      # Set Line style
      histo1.SetLineStyle(0)

      # Set Marker Style
      histo1.SetMarkerStyle(20)
      ErrGraph1.SetMarkerStyle(20)
      histo2.SetMarkerStyle(20)
      ErrGraph2.SetMarkerStyle(20)
      histo3.SetMarkerStyle(20)
      ErrGraph3.SetMarkerStyle(20)
      histo4.SetMarkerStyle(20)
      ErrGraph4.SetMarkerStyle(20)

      # Set Line color
      histo1.SetLineColor(color[0])
      ErrGraph1.SetLineColor(color[0])
      histo2.SetLineColor(color[1])
      ErrGraph2.SetLineColor(color[1])
      histo3.SetLineColor(color[2])
      ErrGraph3.SetLineColor(color[2])
      histo4.SetLineColor(color[3])
      ErrGraph4.SetLineColor(color[3])

      # Set marker color
      histo1.SetMarkerColor(color[0])
      ErrGraph1.SetMarkerColor(color[0])
      histo2.SetMarkerColor(color[1])
      ErrGraph2.SetMarkerColor(color[1])
      histo3.SetMarkerColor(color[2])
      ErrGraph3.SetMarkerColor(color[2])
      histo4.SetMarkerColor(color[3])
      ErrGraph4.SetMarkerColor(color[3])

      # Set range
      histo1.GetYaxis().SetRangeUser(0.0, 1.1)
      histo1.GetXaxis().SetRangeUser(20.0, 300.0)
      if (closer):
        if (multiHistoName == "L1_J75.31ETA49-HLT_j60_320eta490"): histo1.GetXaxis().SetRangeUser(20.0, 250.0)


      # Add legend entry
      legend.AddEntry(histo1,sampleName1, "p")
      legend.AddEntry(histo2,sampleName2, "p")
      legend.AddEntry(histo3,sampleName3, "p")
      legend.AddEntry(histo4,sampleName4, "p")

      # check the histoCounter and therefore if Draw need to be with the option "Same" or without
      if (histoCounter == 0):
        histo1.Draw("hist p")
        MultiGr.Add(ErrGraph1)
#      else:
        #histo1.Draw("Same hist p")
        #MultiGr.Add(ErrGraph1)

      #histo2.Draw("Same hist p")
      #MultiGr.Add(ErrGraph2)

      histo3.Draw("Same hist p")
      MultiGr.Add(ErrGraph3)

      histo4.Draw("Same hist p")
      MultiGr.Add(ErrGraph4)

      histoCounter +=1

      c.Update()

      # Draw legend and ATLAS style label
      legend.Draw("Same")
      AtlasStyle.ATLAS_LABEL(0.5, 0.55, internal = False, preliminary = True, color=1)
      AtlasStyle.myText(0.5, 0.5, 1, "#sqrt{s} = 13 TeV")


      c.Update()

    # Set Title
    l = TLatex(20,1.11, "#font[12]{" + multiHistoName + "}")
    l.Draw("Same")

    # Draw multiGraph
    MultiGr.Draw("P")

    # Write out to file
    MultiGr.Write()

    c.Update()

  # save pdf and eps files
  ## convert names latex friendly
  triggerPDFNamePrae = re.sub("_","-", multiHistoName)
  triggerPDFName =  re.sub(".31ETA","-31ETA", triggerPDFNamePrae)

  c.SaveAs("makeResultPdf/plots/"+"efficienciesPT-" + triggerPDFName + ".pdf")
  c.SaveAs("makeResultPdf/plots/"+"efficienciesPT-" + triggerPDFName + ".eps")


outputfile.Close()



