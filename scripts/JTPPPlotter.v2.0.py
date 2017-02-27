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

inputDataList = ["hist-Data16_run300800_newCut.root", "hist-Data16_run300800_oldCuts.root"]

# please put here the header string in the turnon plots for each file
inputDataTurnonHeader = ["new Cut", "old Cut"]

# Select ATLAS Labeling
ATLAS_Internal = True
ATLAS_Preliminary = False

####### KINEMATICS #######

# Enabling Kinematic plots
DoKinematics = True

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

CriteriaOfHistoCombination = "inputData" # other options: "inputData", "trigger", "jetType", "praefix", "obs", "etaBins"

# Select ratio plots at the buttom of the canvas
# Only available if EnableHistoCombination = True (need other plots to compare with, right? ;) )
# Note: Ratios are always taken with respect to the first histogram in the combined plot
EnableRatioPlots = True

# Set Legend Settings
LegendHeader = "look at me!"
LegendEntryList = ["trigger jet", "offline jet"]

# Line and Marker Style
KLineStyle = 1
KMarkerStyle = 0

######## TURNONS ########

# Enabling Turnon plots
#DoTurnons = False

# Enable Turnon Plotting
PlotTDT             = False
PlotEmulation       = True
PlotBeforePraescale = False #TODO

# Enable Combined Plots
EnableCombinedTurnons = False

# NOTE: If several datasets are listed in inputDataList, the same
# turnons from the different datasets are plotted together in one plot

if (not EnableCombinedTurnons):

    #select if all turnons from logfile should be plotted
    PlotAllTurnonsFromLogfile = False #TODO

    #otherwise create turnon list
    TurnonList=["HLT_j360-HLT_j260","HLT_j110-HLT_j85","HLT_j25-HLT_j15","HLT_j60-HLT_j25"]

    #set the min and the max x-axis range for the selected turnons respectively
    minXList = [250.0,  20.0,  20.0,  20.0]
    maxXList = [600.0, 160.0,  40.0, 110.0]

    #set legend entries
    # Set Legend Settings
    TurnonLegendHeader = "look at Turnons!"
    TurnonLegendEntryList = ["HLT_j360","HLT_j110","HLT_j25","HLT_j60"]

    #rebin factor (has to be a factor of 1200)
    rebinList = [5, 5, 5, 5]

else:

    #create list of name combined plots
    CombinedTurnonsList=["HLT-high", "HLT-low"]

    #set the min and the max x-axis range for the selected turnons respectively
    minXList = [ 20.0,  20.0]
    maxXList = [450.0, 160.0]

    #set legend entries
    TurnonLegendHeader = "look at Turnons!"
    TLegendEntryList = CombinedTurnonsList

    #rebin factor (has to be a factor of 1200)
    rebinList = [5, 5, 5, 5]

    #declare list of contents for every combined plot (do not touch this line)
    ContentCombinedTurnonsList = []

    #select content for combined plot number 1
    ContentCombinedTurnonsList.append = ["HLT_j360-HLT_j260","HLT_j110-HLT_j85"]

    #select content for combined plot number 2
    ContentCombinedTurnonsList.append = ["HLT_j25-HLT_j15","HLT_j60-HLT_j25"]

    #can be continued....

# Line and Marker Style
TLineStyle = 0
TMarkerStyle = [20, 24] # Marker style for different input data

######## GENERAL STYLE ########

# list of colors (feel free to change colors if you are not satisfied)
colorList = [ kRed+2, kBlue+2, kGreen+2, kCyan+2, kMagenta+2, kYellow+2, kRed-2, kBlue+0, kGreen-6, kCyan-2, kMagenta-7, kYellow-7 ]

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
                            histo.SetLineStyle(KLineStyle)
                            histo.SetMarkerStyle(KMarkerStyle)
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
                                        line.Draw("Same"); 
                                        c.Update()


                            # increasing NPlots
                            if (EnableHistoCombination): NPlots += 1


                            if (not EnableHistoCombination):
                                c.SaveAs("makeResultPdf/plots/"+ trigger + "-" + praefix + "-" + jetType + "-" + etaBins + "-" + obs +".pdf")
                                c.SaveAs("makeResultPdf/plots/"+ trigger + "-" + praefix + "-" + jetType + "-" + etaBins + "-" + obs +".eps")

                        if (EnableHistoCombination):
                            # modify name with respect to Criteria
                            subName = trigger + "-" + praefix + "-" + jetType + "-" + etaBins + "-" + obs
                            if (CriteriaOfHistoCombination != "inputData"): 
                                modSubName = subName.replace(var[listSize-1], "")
                            else:
                                modSubName = subName
                            c.SaveAs("makeResultPdf/plots/"+ "Comb" + "-" + modSubName +".pdf")
                            c.SaveAs("makeResultPdf/plots/"+ "Comb" + "-" + modSubName +".eps")


################## CREATING TURNONS  ########################################

# TODO
#ErrGraph1 = TGraphAsymmErrors(histo1, histo1Denum, 'b(1,1) mode')

# Gstyle and ATLAS Style
AtlasStyle.SetAtlasStyle()
gStyle.SetOptStat(0) # get rid of statistics box

# opening all output files
finList = []

for inputData in inputDataList:
    finList.append(TFile.Open(inputData))


# Create PraefixList
EvalTypeList = []

if (PlotTDT): EvalTypeList.append("TDT")
if (PlotEmulation): EvalTypeList.append("Emu")
if (PlotBeforePraescale): EvalTypeList.append("TODO")


for EvalType in EvalTypeList:

    NTurnon = 0
    for TurnonName in TurnonList:

        # create canvas
        c = TCanvas()
        legend = TLegend(0.55,0.50,0.8,0.2)
        legend.SetHeader(TurnonLegendHeader)

        NData = 0
        for fin  in finList:

            # create histo names
            histoName = "TurnOns" + "/" + "effic_pt_" + EvalType + "_" + TurnonName
            histoNameDN = "TurnOns" + "/" + "effic_DENUM_pt_" + EvalType + "_" + TurnonName

            # load histograms from Data
            histo = fin.Get(histoName)
            histoDN = fin.Get(histoNameDN)

            # rebin
            print "NTurnon: " + str(NTurnon)
            histo.Rebin(rebinList[NTurnon])
            histoDN.Rebin(rebinList[NTurnon])

            # Divide to obtain Turnons
            histo.Divide(histo, histoDN, 1.0, 1.0 , "")

            # HISTO COSMETICS

            ### histo style
            histo.SetLineStyle(TLineStyle)
            print "TMarkerStyle[NData]: " + str(TMarkerStyle[NData])
            histo.SetMarkerStyle(TMarkerStyle[NData])
            histo.SetLineColor(colorList[0])
            histo.SetMarkerColor(colorList[0])

            # set range
            histo.GetXaxis().SetRangeUser(minXList[NTurnon], maxXList[NTurnon])
            histo.GetYaxis().SetRangeUser(0.0, 1.2)

            c.Update()

            ### axis style
            histo.GetXaxis().SetTitleOffset(1.4)
            histo.GetYaxis().SetTitleOffset(1.4)

            histo.GetYaxis().SetTitle("entries/width")
            histo.GetXaxis().SetTitle("pt [GeV]")

            c.Update()

            # Drawing
            ## Draw histo
            if (NData == 0): histo.Draw("PE")
            else: histo.Draw("Same PE")

            c.Update()

            ## Draw legend
            legend.AddEntry(histo,TurnonLegendEntryList[NTurnon] + " - " + inputDataTurnonHeader[NData],"p")
            legend.Draw("Same")

            # ATLAS Style Setting
            AtlasStyle.ATLAS_LABEL(0.5, 0.55, internal = ATLAS_Internal, preliminary = ATLAS_Preliminary, color=1)
            AtlasStyle.myText(0.5, 0.5, 1, "#sqrt{s} = 13 TeV")

            NData +=1

            c.Update()

        NTurnon += 1

        # Save Turnons
        c.SaveAs("makeResultPdf/plots/efficienciesPT" + TurnonName +".pdf")
        c.SaveAs("makeResultPdf/plots/efficienciesPT" + TurnonName +".eps")
