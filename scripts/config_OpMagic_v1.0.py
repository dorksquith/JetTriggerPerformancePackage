# ===================================================================================================
# config file for JetTriggerPerformancePackage (JTPP)
#
# Please follow the instructions below in order to set the parameters in the correct way
# and to obtain performance plots of your choice.
# ===================================================================================================

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
# In the following, the initial parameters for JTPP are set.
# Please change ONLY the boolean values (True or False) and the string content
# but NOT the expressions of the initial variables (i.e. the expressions starting with 'm_'
#
c.setalg("InputHandler", { "m_name"                   : "MultijetAlgo",
#
# ============= 0. General Parameters ================================================================================================================
#
# Determine whether the ntuple contains REAL DATA or MC. In case of MC, this option will automatically switch of parameters that exist for real data only
#
                           "m_isData"                                : True,                        
#
# Debug options
#
                           "m_debug"                                 : False,
                           "m_debugInExecute"                        : False,
#
# Determine the name of the offline and trigger branches like "jet_X" and "trigJet_X" (with a capital X!):
#
                           "m_branch_offlineJet"                     : "jet_X",
                           "m_branch_triggerJet"                     : "trigJet_X",
#
# Select how many events you would like to study for this run.
# If '-1' is selected, all events are used.
#
                           "m_doOnlyThisNumberOfEvents"              : 1000,
#
# Specific options for REAL DATA
#
                           "m_checkLArError"                         : True,
#
#
# ============= 1. Kinematics ========================================================================================================================
#
# The following parameters will only affect the kinematic plots
#
# Switch on the generation of Kinematic Plots for offline, trigger and/or truth jets
#
                           "m_doOfflineJetKinematics"                : True,
                           "m_doTriggerJetKinematics"                : True,
                           "m_doTruthJetKinematics"                  : False,
#
# Select triggers which kinematic plots you would like to generate.
# STUDYALL is the option to obtain kinematics for ANY trigger that are included in the ntuple
#
                           "m_TriggerName"                           : "STUDYALL; HLT_j460_a10_lcw_L1HT190-J15.ETA21; HLT_j460_a10_lcw_L1J100; HLT_j25; HLT_j60; HLT_j175; HLT_j360; HLT_j25_320eta490",
#
# Generate leading, subleading and third jet (in pt ordering) and/or the nth jet of each jet trigger
#
                           "m_doLeadSubleadThirdJetKinematics"       : True,
                           "m_doNthJetKinematics"                    : False,
#
# Option to separate kinematic plots in BINS OF ETA
#
                           "m_doKinematicsInBinsOfEta"               : False,
#
# Options to plot kinematics of additional observables:
#
                           "m_doyStar"                               : True,
                           "m_doDeltaPhi"                            : False,
                           "m_doPTBalance"                           : False,
                           "m_doHT"                                  : True,
                           "m_doMHT"                                 : False,
                           "m_doMHTPhi"                              : False,
                           "m_doEMFrac"                              : False,
                           "m_doHECFrac"                             : False,
                           "m_doFracSamplingMax"                     : False,
#
# Specifiy the wanted event cuts for kinematic plots in the string (separate with semicolon, vector in edgy brackets [] and no units).
# For trigger or truth jets, write 'trig' or 'truth' in front of the observable, e.g. trigpt or trig-pt.
# e.g. "|yStar| < 0.6; pt[0] > 150; pt[1] > 50"
#
                           "m_doCuts"                                : True,
                           "m_cutStringKinematics"                   : "trig-pt[0] > 50; |yStar| < 0.2; trig-phi[1] < 2.5",

#
# ============= 2. Response ==========================================================================================================================
#
#
# Switch on the generation of Pt Response Plots for offline, trigger and/or truth jets
#
                           "m_doOfflineTruthResponse"                : False,
                           "m_doTriggerTruthResponse"                : False,
                           "m_doTriggerOfflineResponse"              : False,
#
# Switch on the generation of Mjj Response Plots for offline, trigger and/or truth jets
#
                           "m_doMjjResponseOffVsTruth"               : False,
                           "m_doMjjResponseTrigVsTruth"              : False,
                           "m_doMjjResponseTrigVsOff"                : False,
#
# ============= 3. TurnOns (Trigger Efficiencies) ===================================================================================================
#
#
# Switch on the generation of turnons
#
                           "m_doTurnOns"                             : True,
#
# Select triggers which turnon plots you would like to generate.
# The syntax for turnon plots is 'probe-trigger'-'reference-trigger', e.g. HLT_j25-HLT_j15
#
                           "m_TurnOnName"                            : "HLT_j25/HLT_j15;  HLT_j60/HLT_j25; HLT_j110/HLT_j85; HLT_j360/HLT_j260;  HLT_j380/HLT_j260; HLT_j400/  HLT_j260",
#
# Select the STRATEGY of turnOn generation.
# The available options are the Trigger Decision Tool (TDT), Emulation or Trigger Before Preascale. Note that more than one strategy can be selected
# and will be created independently of each other
#
                           "m_useTriggerDecisionTool"                : True,
                           "m_useEmulation"                          : False,
                           "m_useTriggerBeforePraescale"             : False,
#
# For the turnon event selection, following cuts on offline jets are applied by default:
# HLT:                       L1:
# |eta| < 2.8 (central)       |eta| < 2.6 (central)
# |eta| > 3.6 (fwd)           |eta| > 3.6 (fwd)
# |timing| < 10.0             |timing| < 10.0
#
# Additionally, you can define the event selection for turnons here:
#
                           "m_cutStringTurnons"                      : "|pt[0]| < 50",

# TODO CUTS
                           #cuts:
# CLEANING                           "m_doCleaning"                            : False,
                           } )

