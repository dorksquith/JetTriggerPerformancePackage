#!/bin/python
########################################################################
# The JTPPPlotter by Edgar Kellermann (2016/03/14)
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
fin1 = TFile.Open("Turnons_ChristiansRootFiles/r29926_PhysicsMain00297041.root")
fin2 = TFile.Open("hist-Data16_AOD_run297041/hist-Data16_AOD_run297041.root")

## sample name in the legend
sampleName1 = "Christian"
sampleName2 = "Edgar"

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
multiHistoList = ["HLT-singleJetTriggers", "L1-singleJetTriggers"]

# -------------------------- TurnOns ----------------------------------

outputfile.cd()
c = TCanvas("c", sampleName1, 1)

# loop over all multi histograms
for multiHistoName in multiHistoList:

  ## Legend
  legend = TLegend(0.5,0.40,0.75,0.25)
  legend.SetFillColor(0)
  legend.SetTextFont(1)
  legend.SetHeader("emulated")
  legend.SetFillStyle(0)


  # multigraph for errors
  MultiGr = TMultiGraph()


  histoCounter = 0

  # Loop over all turnons
  for turnOn in turnOnList:


    print "Filling... " + multiHistoName
    if (multiHistoName == "HLT-singleJetTriggers"):
      if (turnOn == "HLT_j25-HLT_j15"): print "using... " + turnOn
      elif (turnOn == "HLT_j60-HLT_j25"): print "using... " + turnOn
      elif (turnOn == "HLT_j110-HLT_j85"): print "using... " + turnOn
      else: continue

    if (multiHistoName == "L1-singleJetTriggers"):
      if (turnOn == "L1_J15-HLT_j15"): print "using... " + turnOn
#      elif (turnOn == "L1_J40-HLT_j45"): print "using... " + turnOn
      elif (turnOn == "L1_J40-HLT_j25"): print "using... " + turnOn
      else: continue

    # Loopt over EvaluationType
    for EvalType in ["Emu"]: # ["Emu", "TDT"]:

      # TEfficiency name in Christian's file
      effic1Name1 = turnOn.rsplit('-',1)[0].replace("J", "j")
      print "effic1Name1: " + effic1Name1

      # histo name in Edgar's files
      histoName2 = "TurnOns" + "/" + "effic_pt_" + EvalType + "_" + turnOn 
      histoName2Denum = "TurnOns" + "/" + "effic_DENUM_pt_" + EvalType + "_" + turnOn

      # load TEfficiencies
      effic1 = fin1.Get(effic1Name1)

      # set correct statistical bars according to Christian
      effic1.SetStatisticOption(TEfficiency.kBUniform);
      effic1.SetPosteriorMode();

      # load histos
      histo2 = fin2.Get(histoName2)
      histo2Denum = fin2.Get(histoName2Denum)

      # create ErrGraph
      ErrGraph2 = TGraphAsymmErrors(histo2, histo2Denum, 'b(1,1) mode')

      # Divide histograms
      histo2.Divide(histo2, histo2Denum, 1.0, 1.0 , "")


      # get rid of non fully efficient region
      if (turnOn == "HLT_j110-HLT_j85"):
        for n in range (0, 37):
          histo2.SetBinContent(n, -1)
          ErrGraph2.SetPoint(n, 1000, 1000)

      if (turnOn == "HLT_j60-HLT_j25"):
        for n in range (0, 13):
          histo2.SetBinContent(n, -1)
          ErrGraph2.SetPoint(n, 1000, 1000)

      # get rid of non fully efficient region in Christian's plots
      if (effic1Name1 == "HLT_j110"):
        for n in range (0, 14):
          effic1.SetPassedEvents(n, 0)

      if (effic1Name1 == "HLT_j60"):
        for n in range (0, 12):
          effic1.SetPassedEvents(n, 0)

      if (effic1Name1 == "L1_j15"):
        for n in range (0, 4):
          effic1.SetPassedEvents(n, 0)

      if (effic1Name1 == "L1_j40"):
        for n in range (0, 11):
          effic1.SetPassedEvents(n, 0)

      # Set Axis Titles
      histo2.GetXaxis().SetTitle("offline pt [GeV]")
      histo2.GetYaxis().SetTitle("efficiency")

      # Set Title offset
      histo2.GetXaxis().SetTitleOffset(1.4)
      histo2.GetYaxis().SetTitleOffset(1.4)

      # Set Line style
      effic1.SetLineStyle(0)
      histo2.SetLineStyle(0)

      # Set Marker Style
      effic1.SetMarkerStyle(20)
      histo2.SetMarkerStyle(24)
      ErrGraph2.SetMarkerStyle(24)

      # Set Line color
      effic1.SetLineColor(color[histoCounter])
      histo2.SetLineColor(color[histoCounter])
      ErrGraph2.SetLineColor(color[histoCounter])

      # Set marker color
      effic1.SetMarkerColor(color[histoCounter])
      histo2.SetMarkerColor(color[histoCounter])
      ErrGraph2.SetMarkerColor(color[histoCounter])

      # Set range
      histo2.GetYaxis().SetRangeUser(0.0, 1.1)
      histo2.GetXaxis().SetRangeUser(20.0, 600.0)
      if (closer):
        if (multiHistoName == "HLT-singleJetTriggers"): histo2.GetXaxis().SetRangeUser(0.0, 250.0)
        if (multiHistoName == "L1-singleJetTriggers"): histo2.GetXaxis().SetRangeUser(0.0, 300.0)


      # Add legend entry
      legend.AddEntry(effic1,turnOn.rsplit('-',1)[0]  + " " + sampleName1, "p")
      legend.AddEntry(histo2,turnOn.rsplit('-',1)[0]  + " " + sampleName2, "p")

      # check the histoCounter and therefore if Draw need to be with the option "Same" or without
      if (histoCounter == 0):
        histo2.Draw("hist p")
        MultiGr.Add(ErrGraph2)
      else:
        histo2.Draw("Same hist p")
        MultiGr.Add(ErrGraph2)

      #MultiGr.Add(effic1)
      effic1.Draw("Same p")

      histoCounter +=1

      c.Update()

      # Draw legend and ATLAS style label
      legend.Draw("Same")
      AtlasStyle.ATLAS_LABEL(0.5, 0.55, internal = False, preliminary = True, color=1)
      AtlasStyle.myText(0.5, 0.5, 1, "#sqrt{s} = 13 TeV")


      c.Update()

    # Set Title
    l = TLatex(20,1.11, "#font[12]{" + sampleName1 + "}")
    l.Draw("Same")

    # Draw multiGraph
    MultiGr.Draw("P")

    # Write out to file
    MultiGr.Write()

    c.Update()

# save pdf and eps files
## convert names latex friendly
#triggerPDFNamePrae = re.sub("_","-", turnOn)
#triggerPDFName =  re.sub(".31ETA","-31ETA", triggerPDFNamePrae)

  c.SaveAs("makeResultPdf/plots/"+"efficienciesPT-" + multiHistoName + ".pdf")
  c.SaveAs("makeResultPdf/plots/"+"efficienciesPT-" + multiHistoName + ".pdf")


outputfile.Close()



