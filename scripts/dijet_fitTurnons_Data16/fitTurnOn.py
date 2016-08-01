#!/usr/bin/env python

from ROOT import *
import math
from scipy.optimize import fsolve
import AtlasStyle

AtlasStyle.SetAtlasStyle()

def find_turn_on(gname,satpoint,xmin,xmax):

    f = TFile.Open("JetTriggerPerformancePlots.root")
    g = f.Get(gname)

    g.Draw("ap")

    # nadia function

    # func = TF1.TF1("sigm","50*TMath::Erf( (x-[1])/sqrt(2*[2]*[2]) )*(1-[0])+50*(1-[0])",xmin,xmax)
    # func.SetParameter(0,0) # overall normalization
    # func.SetParameter(1,(xmax+xmin)/2.) # position of half-max
    # func.SetParameter(2,0.001) # slope of turn-on. 100 still ok
    # func.SetParLimits(0,0,1)
    # func.SetParLimits(1,xmin,xmax)
    # func.SetParLimits(2,0.000001,1000)
    # g.Fit(func,"","",xmin,xmax)

    # sigmoid function

    # func = TF1.TF1("sigm","[0]/(1+exp(-2*(x-[1])/[2]))",xmin,xmax)
    # func.SetParameter(0,100) # overall normalization
    # func.SetParameter(1,(xmax+xmin)/2.) # position of half-max
    # func.SetParameter(2,0.001) # slope of turn-on. 100 still ok
    # func.SetParLimits(0,0,100)
    # func.SetParLimits(1,xmin,xmax)
    # func.SetParLimits(2,0.000001,1000)
    # g.Fit(func,"","",xmin,xmax)

    # double sigmoid

    func = TF1.TF1("sigm","([0]*(1 - 1/(1 + exp(([1] + [2]/2. - x)/[4]))))/(1 + exp(([1] - [2]/2. - x)/[3]))",xmin,xmax)
    func.SetParameter(0,200) # overall normalization
    func.SetParameter(1,(xmax+xmin)/3.) # position of half-max
    func.SetParameter(2,0) # position of second
    func.SetParameter(3,60)  # slope of turn-on 1.
    func.SetParameter(4,1000)  # slope of turn-on 2.
    func.SetLineColor(kRed)
    g.Fit(func,"","",xmin,xmax)

    # find saturation point (e.g. satpoint=99.5)
    fitresult = lambda x: (func.Eval(x)-satpoint)
    x_init = xmax
    x_satpoint = fsolve(fitresult,x_init)
    y_satpoint = fitresult(x_satpoint)+satpoint
    print "x_satpoint: " + str(x_satpoint)
    print "y_satpoint: " + str(y_satpoint)

    line = TLine()
    text = TText()
    line.SetLineStyle(kDashed)
    line.SetLineColor(kRed)
    line.DrawLine(x_satpoint,0,x_satpoint,100)
    text.SetTextColor(kRed)
    text.SetTextFont(1)
    text.SetTextSize(0.04)
    text.DrawText(x_satpoint*0.9,50,"y:%6.4f"%(y_satpoint))
    text.DrawText(x_satpoint*0.9,45,"x:%6.4f"%(x_satpoint))

    g.GetXaxis().SetRangeUser(xmin*0.6,xmax*1.1)
    gPad.Print( (gname+"_goal%04d_x%05.1f_y%05.1f.pdf") % (satpoint*1,x_satpoint[0],y_satpoint) )

    return x_satpoint[0]

    #func.SetParLimits(1,xmin,xmax)
    #func.SetParLimits(2,0.000001,1000)


hFile = TFile.Open("JetTriggerPerformancePlots.root")
histKeys = hFile.GetListOfKeys()
for hkey in histKeys:
    hname = hkey.GetName()
    if "HLT_j25-" in hname:
        find_turn_on( hname , 99.5, 20, 80 )
    if "HLT_j60-" in hname:
        find_turn_on( hname , 99.5, 40, 100 )
    if "HLT_j110-" in hname:
      find_turn_on( hname , 99.5, 90, 150 )
    if "HLT_j360-" in hname:
      find_turn_on( hname , 99.5, 250, 550 )
    if "HLT_j380-" in hname:
        find_turn_on( hname , 99.5, 300, 500 )
    if "HLT_j400-" in hname:
        find_turn_on( hname , 99.5, 300, 500 )


