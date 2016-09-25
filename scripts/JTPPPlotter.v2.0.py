###########################################################################
#
# This is JTPPPlotter.v2.0.py providing plot mechanisms for JTPP.
#
# by Edgar Kellermann (edgar.kellermann@cern.ch)
###########################################################################

from ROOT import *
import math
import sys
import re
from array import array
import AtlasStyle

################### CONFIGURATIONS ########################################

# Select input data:

#inputDataList = ["hist-Data16_run300800.root", "hist-Data16_run300800_oldCuts.root"]
inputDataList = ["hist-Data16_run300800_oldCuts.root"]

# Select ATLAS Labeling
ATLAS_Internal = True
ATLAS_Preliminary = False

####### KINEMATICS #######

# Enabling Kinematic plots
DoKinematics = False

# Select wanted plots:

## Select plots from logfile (i.e. all provided plots)
ReadLogfileForKinematics = False

if (ReadLogfileForKinematics):
    print "Reading Logfile...."
    # TODO Implement reading logfile

## Otherwise pick wanted plots and but them into the lists
else:

    # jetType (jet = offline/reco jets)
    jetTypeList = ["trigJet", "jet"]
    #["jet", "trigJet", "truthJet"]

    # all jets, leading, sublead, third, nth jets
    praefixList = [""]
    # ["","lead", "sublead", "third", "nth"]

    # observables
    obsList = ["E","phi"]
    # ["E", "pt", "phi", "eta", "mjj", "yStar", "mjj", "m23", "yStar", "deltaPhi", "pTBalance", "MHT", "MHTPhi", "EMFrac", "HECFrac", "FracSamplingMax"]

    # etaBins # TODO Implement proper list
    etaBinsList = [""]

    # triggers of interest
    triggerList = ["HLT_j360"]

# Select if histos should be plotted together on one canvas and  the criteria for combination
EnableHistoCombination = True

CriteriaOfHistoCombination = "jetType" # other options: "inputData", "trigger", "jetType", "praefix", "obs", "etaBins"

# Select ratio plots at the buttom of the canvas
# Only available if EnableHistoCombination = True (need other plots to compare with, right? ;) )
# Note: Ratios are always taken with respect to the first histogram in the combined plot
EnableRatioPlots = True

# Set Legend Settings
LegendHeader = "look at me!"
LegendEntryList = ["trigger jet", "offline jet"]

######## TURNONS ########

# Enable Turnon Plotting
PlotTDT             = False
PlotEmulation       = True
PlotBeforePraescale = False #TODO

# Enable Combined Plots
EnableCombinedTurnons = False

if (not EnableCombinedTurnons):

    #select if all turnons from logfile should be plotted
    PlotAllTurnonsFromLogfile = False #TODO

    #otherwise create turnon list
    TurnonList=["HLT_j360-HLT_j260","HLT_j110-HLT_j85","HLT_j25-HLT_j15","HLT_j60-HLT_j25"]


else:

    #create list of name combined plots
    CombinedTurnonsList=["HLT-high", "HLT-low"]

    #declare list of contents for every combined plot (do not touch this line)
    ContentCombinedTurnonsList = []

    #select content for combined plot number 1
    ContentCombinedTurnonsList.append = ["HLT_j360-HLT_j260","HLT_j110-HLT_j85"]

    #select content for combined plot number 2
    ContentCombinedTurnonsList.append = ["HLT_j25-HLT_j15","HLT_j60-HLT_j25"]

    #can be continued....

######## STYLE ########

# list of colors (feel free to change colors if you are not satisfied)
colorList = [ kRed+2, kBlue+2, kGreen+2, kCyan+2, kMagenta+2, kYellow+2, kRed-2, kBlue+0, kGreen-6, kCyan-2, kMagenta-7, kYellow-7 ]

# Line and Marker Style
LineStyle = 1
MarkerStyle = 0

###############################################################################

################## CREATING KINEMATICS ########################################

if (DoKinematics):

    # Initial style settings
    AtlasStyle.SetAtlasStyle()
    gStyle.SetOptStat(0) # get rid of statistics box

    # opening all output files
    finList = []

    for inputData in inputDataList:
        finList.append(TFile.Open(inputData))

    # Looping in a proper way
    listSize = 6

    list = [inputDataList] *listSize
    var = ["ThereWillBeString"] * listSize

    ## permuting lists depending on the combination criteria (since the criteria has to be in the last loop)
    i = -1
    if (CriteriaOfHistoCombination == "etaBins"):   i = 0
    if (CriteriaOfHistoCombination == "obs"):       i = 1
    if (CriteriaOfHistoCombination == "praefix"):   i = 2
    if (CriteriaOfHistoCombination == "jetType"):   i = 3
    if (CriteriaOfHistoCombination == "trigger"):   i = 4
    if (CriteriaOfHistoCombination == "inputData"): i = 5
    if (not EnableHistoCombination):                i = 0

    # assert
    if (i == -1): print "ERROR: Combination Criteria cannot be read"; # TODO ASSERT

    list[(0+i)%listSize] = finList
    list[(1+i)%listSize] = triggerList
    list[(2+i)%listSize] = jetTypeList
    list[(3+i)%listSize] = praefixList
    list[(4+i)%listSize] = obsList
    list[(5+i)%listSize] = etaBinsList

    for var[0] in list[0]:
        for var[1] in list[1]:
            for var[2] in list[2]:
                for var[3] in list[3]:
                    for var[4] in list[4]:

                        #create Canvas
                        if (EnableHistoCombination):
                            c = TCanvas()
                            NPlots = 0
                            legend = TLegend(0.5,0.60,0.75,0.75)
                            legend.SetHeader(LegendHeader)

                            if (EnableRatioPlots):
                                #making the pads for the histoRatio:
                                pad1 = TPad("pad1","pad1",0.05,0.30,1,1);
                                pad2 = TPad("pad2","pad2",0.05,0.05,1,0.30);
                                pad1.SetTopMargin(0.02);
                                pad2.SetTopMargin(0.0);
                                pad1.SetBottomMargin(0.0);
                                pad2.SetBottomMargin(0.30);

                        for var[5] in list[5]:

                            #create Canvas
                            if (not EnableHistoCombination):
                                c = TCanvas()
                                NPlots = 0
                                legend = TLegend(0.5,0.60,0.75,0.75)
                                legend.SetHeader(LegendHeader)

                            # identify what the var[] list actually contains:
                            fin = var[(0+i)%listSize]
                            trigger   = var[(1+i)%listSize]
                            jetType   = var[(2+i)%listSize]
                            praefix   = var[(3+i)%listSize]
                            obs       = var[(4+i)%listSize]
                            etaBins   = var[(5+i)%listSize]

                            # no lead, sublead, third for scalar observables:
                            if (((obs == "mjj")|(obs == "m23")|(obs == "yStar")|(obs == "deltaPhi")|(obs == "pTBalance")|(obs == "MHT")|(obs == "MHTPhi")) & (praefix != "")): continue

                            # loading histogram
                            histoName = trigger + "/" + jetType + "Kinematics/" + praefix + jetType + etaBins +  "_" + obs

                            # opening file
                            #if (NPlots == 0): fin = TFile.Open(inputData)

                            # loading histo
                            histo = fin.Get(histoName)

                            # HISTO COSMETICS

                            ### histo style
                            histo.SetLineStyle(LineStyle)
                            histo.SetMarkerStyle(MarkerStyle)
                            histo.SetLineColor(colorList[NPlots])
                            histo.SetMarkerColor(colorList[NPlots])

                            ### axis style
                            histo.GetXaxis().SetTitleOffset(1.4)
                            histo.GetYaxis().SetTitleOffset(1.4)

                            histo.GetYaxis().SetTitle("entries/width")
                            if ((obs == "pt")| (obs == "E")| (obs == "mjj")| (obs == "m23")| (obs == "MHT")):
                                histo.GetXaxis().SetTitle(obs + " [GeV]")

                            c.Update()

                            # draw pads and change to pad1
                            if (EnableHistoCombination & EnableRatioPlots & (NPlots == 0)):
                                pad1.Draw()
                                pad2.Draw()
                                pad1.cd()

                            # drawing plots
                            if (NPlots == 0): histo.Draw("HIST")
                            else: histo.Draw("Same HIST")

                            c.Update()

                            # legend
                            legend.AddEntry(histo,LegendEntryList[NPlots],"l")
                            legend.Draw("Same")

                            # ATLAS Style Setting
                            AtlasStyle.ATLAS_LABEL(0.5, 0.85, internal = ATLAS_Internal, preliminary = ATLAS_Preliminary, color=1)
                            AtlasStyle.myText(0.5, 0.8, 1, "#sqrt{s} = 13 TeV")

                            # create ratio plots in pad2
                            if (EnableHistoCombination & EnableRatioPlots):

                                ## save reference plot, i.e. NPlots == 0
                                if (NPlots == 0): refHisto = histo

                                ## else plot ratio histos
                                else:
                                    # change to pad2
                                    pad2.cd()
                                    pad2.SetLogy(False)

                                    ## create and setup histoRatio
                                    histoRatio=histo.Clone()
                                    histoRatio.Sumw2()
                                    histoRatio.Divide(refHisto)

                                    histoRatio.GetXaxis().SetTitleSize(0.07)
                                    histoRatio.GetXaxis().SetLabelFont(42)
                                    histoRatio.GetXaxis().SetLabelSize(0.1)
                                    histoRatio.GetXaxis().SetTitleSize(0.1)
                                    histoRatio.GetYaxis().SetRangeUser(0.5,1.5)
                                    histoRatio.GetYaxis().SetTitle("Ratio")
                                    histoRatio.GetYaxis().SetTitleOffset(0.7)
                                    histoRatio.GetYaxis().SetLabelFont(42)
                                    histoRatio.GetYaxis().SetLabelSize(0.1)
                                    histoRatio.GetYaxis().SetTitleSize(0.08)

                                    if (NPlots == 1): histoRatio.Draw("hist")
                                    else: histoRatio.Draw("Same hist")
                                    c.Update()

                                    # black lines at 1.0 in ratio plot
                                    for yLine in [0.9, 1.0, 1.1]:
                                        line = TLine(histoRatio.GetXaxis().GetXmin(), yLine, histoRatio.GetXaxis().GetXmax(), yLine);
                                        line.SetLineColor(kBlack);
                                        line.SetLineStyle(3);
                                        #line.Draw("Same"); WARNING
                                        c.Update()




                            # increasing NPlots
                            if (EnableHistoCombination): NPlots += 1


                            if (not EnableHistoCombination):
                                c.SaveAs("makeResultPdf/plots/"+ trigger + "-" + praefix + "-" + jetType + "-" + etaBins + "-" + obs +".pdf")
                                c.SaveAs("makeResultPdf/plots/"+ trigger + "-" + praefix + "-" + jetType + "-" + etaBins + "-" + obs +".eps")

                        if (EnableHistoCombination):
                            # modify name with respect to Criteria
                            subName = trigger + "-" + praefix + "-" + jetType + "-" + etaBins + "-" + obs
                            modSubName = subName.replace(var[listSize-1], "")
                            c.SaveAs("makeResultPdf/plots/"+ "Comb" + "-" + modSubName +".pdf")
                            c.SaveAs("makeResultPdf/plots/"+ "Comb" + "-" + modSubName +".eps")


################## CREATING TURNONS  ########################################



