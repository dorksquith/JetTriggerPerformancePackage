import ROOT
from xAH_config import xAH_config

c = xAH_config()

GRL = "$ROOTCOREBIN/data/ZprimeDM/data15_13TeV.periodAllYear_DetStatus-v64-pro19_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml"

#
#  Data Config
#
if not args.is_MC:
    applyGRL           = True

#
#  MC Config
#
else:
    applyGRL           = False


doTruthOnly = True

#
# Process Ntuple
#
c.setalg("InputHandler", { "m_name"                   : "MultijetAlgo",
                           #triggers of Interest
                           "m_TriggerName"                           : "STUDYALL; HLT_j25; HLT_j60; HLT_j175; HLT_j360; HLT_j25_320eta490; HLT_j60_320eta490; HLT_j110_320eta490; HLT_j260_320eta490; HLT_4j25; HLT_4j45; HLT_5j25; HLT_5j45; HLT_6j25; HLT_7j25; L1_J75.31ETA49",
#                           "m_TriggerName"                           : "STUDYALL; HLT_j25",
                           #debug
                           "m_debug"                                 : True,
                           "m_debugInExecute"                        : False,
                           #Data or MC:
                           "m_isData"                                : True,
                           #For Data samples:
                           "m_checkLArError"                         : True,
                           #Kinematics:
                           "m_doOfflineJetKinematics"                : False,
                           "m_doTriggerJetKinematics"                : False,
                           "m_doTruthJetKinematics"                  : False,
                           #More setups for kinematic
                           "m_doLeadSubleadThirdJetKinematics"       : False,
                           "m_doNthJetKinematics"                    : False,
                           "m_doKinematicsInBinsOfEta"               : False,
                           "m_doEtainBinsOfPt"                       : False,
                           #More Kinematics
                           "m_doyStar"                               : False,
                           "m_doDeltaPhi"                            : False,
                           "m_doPTBalance"                           : False,
                           "m_doMHT"                                 : False,
                           "m_doMHTPhi"                              : False,
                           "m_doEMFrac"                              : False,
                           "m_doHECFrac"                             : False,
                           "m_doFracSamplingMax"                     : False,
                           #ptResponse:
                           "m_doOfflineTruthResponse"                : False,
                           "m_doTriggerTruthResponse"                : False,
                           "m_doTriggerOfflineResponse"              : False,
                           #mjjResponse:
                           "m_doMjjResponseOffVsTruth"               : False,
                           "m_doMjjResponseTrigVsTruth"              : False,
                           "m_doMjjResponseTrigVsOff"                : False,
                           #Turn-ons
                           "m_doTurnOns"                             : True,
                           "m_TurnOnName"                            : "HLT_j25-HLT_j15; HLT_j60-HLT_j25; HLT_j110-HLT_j85; HLT_j360-HLT_j260; HLT_j380-HLT_j260; HLT_j400-HLT_j260; L1_J15-HLT_j15; L1_J20-HLT_j25; L1_J40-HLT_j25; L1_J50-HLT_j25; L1_J40-HLT_j45; L1_J100-HLT_j85; L1_J75-HLT_j45; HLT_j110-HLT_j45; HLT_3j175-HLT_j110; HLT_4j85-HLT_j85; HLT_6j45-HLT_j60; HLT_6j45-HLT_j45; L1_3J40-HLT_j60; L1_4J15-HLT_j45; HLT_j25_320eta490-HLT_j15_320eta490; HLT_j60_320eta490-HLT_j25_320eta490; L1_J15.31ETA49-HLT_j15_320eta490; L1_J75.31ETA49-HLT_j60_320eta490; L1_J75.31ETA49-HLT_j45_320eta490; L1_J20.31ETA49-HLT_j15_320eta490; L1_J50.31ETA49-HLT_j45_320eta490",
                           "m_useTriggerDecisionTool"                : False,
                           "m_useEmulation"                          : True,
                           "m_useTriggerBeforePraescale"             : False,
                           #cuts:
                           "m_doCuts"                                : False,
                           "m_doCleaning"                            : False,
                           "m_doOnlyThisNumberOfEvents"              : -1,
                           } )

