# ===================================================================================================
# config file for JetTriggerPerformancePackage (JTPP)
#
# Please follow the instructions below in order to set the parameters in the correct way
# and to obtain performance plots of your choice.
# ===================================================================================================

import ROOT
from xAH_config import xAH_config

c = xAH_config()


#
#  Data Config
#
if not args.is_MC:
    applyGRL           = True
    isData             = True

#
#  MC Config
#
else:
    applyGRL           = False
    isData             = False

#
# In the following, the initial parameters for JTPP are set.
# Please change ONLY the boolean values (True or False) and the string content
# but NOT the expressions of the initial variables (i.e. the expressions starting with 'm_'
#
c.setalg("InputHandler", { "m_name"                   : "MultijetAlgo",
#
# ============= 0. General Parameters ================================================================================================================
#
# Determine whether the ntuple contains REAL DATA or MC. In case of MC, this option will automatically switch off parameters that exist for real data only.
# By Default, this option is set by the tag --isMC when xAH_run.py is used.
#
                           "m_isData"                                : isData,
#
# GRL options
#
                           "m_applyGRL"                              : True,
                           "m_GRLxml"                                : "data16_13TeV.periodAllYear_DetStatus-v83-pro20-15_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml",
#
# Debug options
#
                           "m_debug"                                 : False,
                           "m_debugInExecute"                        : False,
                           "m_coutPassedTriggers"                    : False,
#
# Determine the name of the offline and trigger branches like "jet_X" and "trigJet_X" (with a capital X!):
#
                           "m_branch_offlineJet"                     : "jet_X",
                           "m_branch_triggerJet"                     : "trigJet_X",
#
# Select triggers which KINEMATIC AND RESPONSE plots you would like to generate.
# STUDYALL is the option to obtain kinematics for ANY trigger that are included in the ntuple.
#
                           "m_TriggerName"                           : "STUDYALL; HLT_j460_a10_lcw_L1HT190-J15.ETA21; HLT_j460_a10_lcw_L1J100; HLT_j25; HLT_j60; HLT_j175; HLT_j360; HLT_j25_320eta490",

#
# Specifiy the wanted event cuts for ALL plots in the string (separate with semicolon, vector in edgy brackets [] and no units).
# For trigger or truth jets, write 'trig' or 'truth' in front of the observable, e.g. trigpt or trig-pt.
# e.g. "|yStar| < 0.6; pt[0] > 150; pt[1] > 50"
#
                           "m_doCuts"                                : False,
                           "m_cutStringKinematics"                   : "trig-pt[0] > 50; |yStar| < 0.2; trig-phi[1] < 2.5",
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
# Enable the generation of Kinematic Plots for offline, trigger and/or truth jets
#
                           "m_doOfflineJetKinematics"                : True,
                           "m_doTriggerJetKinematics"                : True,
                           "m_doTruthJetKinematics"                  : False,
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
                           "m_doMjj"                                 : True,
                           "m_doM23"                                 : True,
                           "m_doyStar"                               : True,
                           "m_doDeltaPhi"                            : False,
                           "m_doPTBalance"                           : False,
                           "m_doHT"                                  : True, #is always calculated by JTPP
                           "m_doMHT"                                 : False,
                           "m_doMHTPhi"                              : False,
                           "m_doEMFrac"                              : False,
                           "m_doHECFrac"                             : False,
                           "m_doFracSamplingMax"                     : False,
#
# Option to calculate Mjj, M23 and yStar on its own in the whole code instead of reading it from the ntuple.
#
                           "m_calculateMjj"                          : True,
                           "m_calculateM23"                          : True,
                           "m_calculateYStar"                        : True,
#
# Binning for the kinematic plots:
# elements have to be separated by kommas ","
#
                           "m_Kinematic_E_Binning"                   : "15. ,20. ,25. ,35. ,45. ,55. ,70. ,85. ,100. ,116. ,134. ,152. ,172. ,194. ,216. ,240. ,264. ,290. ,318. ,346.,376.,408.,442.,478.,516.,556.,598.,642.,688.,736.,786.,838.,894.,952.,1012.,1076.,1162.,1250.,1310.,1420.,1530.,1750.,1992.,2242.,2500.,2850.,3200.,3600.,4000.,4600.",
                           "m_Kinematic_pt_Binning"                  : "15. ,20. ,25. ,35. ,45. ,55. ,70. ,85. ,100. ,116. ,134. ,152. ,172. ,194. ,216. ,240. ,264. ,290. ,318. ,346.,376.,408.,442.,478.,516.,556.,598.,642.,688.,736.,786.,838.,894.,952.,1012.,1076.,1162.,1250.,1310.,1420.,1530.,1750.,1992.,2242.,2500.,2850.,3200.,3600.,4000.,4600.",
                           "m_Kinematic_phi_Binning"                 : "-3.5, -3.36, -3.22, -3.08, -2.94, -2.8, -2.66, -2.52, -2.38, -2.24, -2.1, -1.96, -1.82, -1.68, -1.54, -1.4, -1.26, -1.12, -0.98, -0.84, -0.7, -0.56, -0.42, -0.28, -0.14, 0.0, 0.14, 0.28, 0.42, 0.56, 0.7, 0. 84, 0.98, 1.12, 1.26, 1.4, 1.54, 1.68, 1.82, 1.96, 2.1, 2.24, 2.38, 2.52, 2.66, 2.8, 2.94, 3.08, 3.22, 3.36, 3.5",
                           "m_Kinematic_eta_Binning"                 : "-4.9,-3.1,-2.8,-2.1,-1.8,-1.5,-1.2,-0.8,0.0,0.8,1.2,1.5,1.8,2.1,2.8,3.1,4.9",
#
                           "m_Kinematic_mjj_Binning"                 : "50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 203.0, 216.0, 229.0, 243.0, 257.0, 272.0, 287.0, 303.0, 319.0, 335.0, 352.0, 369.0, 387.0, 405.0, 424.0, 443.0, 462.0, 482.0, 502.0, 523.0, 544.0, 566.0, 588.0, 611.0, 634.0, 657.0, 681.0, 705.0, 730.0, 755.0, 781.0, 807.0, 834.0, 861.0, 889.0, 917.0, 946.0, 976.0, 1006.0, 1037.0, 1068.0, 1100.0, 1133.0, 1166.0, 1200.0, 1234.0, 1269.0, 1305.0, 1341.0, 1378.0, 1416.0, 1454.0, 1493.0, 1533.0, 1573.0, 1614.0, 1656.0, 1698.0, 1741.0, 1785.0, 1830.0, 1875.0, 1921.0, 1968.0, 2016.0, 2065.0, 2114.0, 2164.0, 2215.0, 2267.0, 2320.0, 2374.0, 2429.0, 2485.0, 2542.0, 2600.0, 2659.0, 2719.0, 2780.0, 2842.0, 2905.0, 2969.0, 3034.0, 3100.0, 3167.0, 3235.0, 3305.0, 3376.0, 3448.0, 3521.0, 3596.0, 3672.0, 3749.0, 3827.0, 3907.0, 3988.0, 4070.0, 4154.0, 4239.0, 4326.0, 4414.0, 4504.0, 4595.0, 4688.0, 4782.0, 4878.0, 4975.0, 5074.0, 5175.0, 5277.0, 5381.0, 5487.0, 5595.0, 5705.0, 5817.0, 5931.0, 6047.0, 6165.0, 6285.0, 6407.0, 6531.0, 6658.0, 6787.0, 6918.0, 7052.0, 7188.0, 7326.0, 7467.0, 7610.0, 7756.0, 7904.0, 8055.0, 8208.0, 8364.0, 8523.0, 8685.0, 8850.0, 9019.0, 9191.0, 9366.0, 9544.0, 9726.0, 9911.0, 10100.0, 10292.0, 10488.0, 10688.0, 10892.0, 11100.0, 11312.0, 11528.0, 11748.0, 11972.0, 12200.0, 12432.0, 12669.0, 12910.0, 13156.0",
                           "m_Kinematic_m23_Binning"                 : "50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 203.0, 216.0, 229.0, 243.0, 257.0, 272.0, 287.0, 303.0, 319.0, 335.0, 352.0, 369.0, 387.0, 405.0, 424.0, 443.0, 462.0, 482.0, 502.0, 523.0, 544.0, 566.0, 588.0, 611.0, 634.0, 657.0, 681.0, 705.0, 730.0, 755.0, 781.0, 807.0, 834.0, 861.0, 889.0, 917.0, 946.0, 976.0, 1006.0, 1037.0, 1068.0, 1100.0, 1133.0, 1166.0, 1200.0, 1234.0, 1269.0, 1305.0, 1341.0, 1378.0, 1416.0, 1454.0, 1493.0, 1533.0, 1573.0, 1614.0, 1656.0, 1698.0, 1741.0, 1785.0, 1830.0, 1875.0, 1921.0, 1968.0, 2016.0, 2065.0, 2114.0, 2164.0, 2215.0, 2267.0, 2320.0, 2374.0, 2429.0, 2485.0, 2542.0, 2600.0, 2659.0, 2719.0, 2780.0, 2842.0, 2905.0, 2969.0, 3034.0, 3100.0, 3167.0, 3235.0, 3305.0, 3376.0, 3448.0, 3521.0, 3596.0, 3672.0, 3749.0, 3827.0, 3907.0, 3988.0, 4070.0, 4154.0, 4239.0, 4326.0, 4414.0, 4504.0, 4595.0, 4688.0, 4782.0, 4878.0, 4975.0, 5074.0, 5175.0, 5277.0, 5381.0, 5487.0, 5595.0, 5705.0, 5817.0, 5931.0, 6047.0, 6165.0, 6285.0, 6407.0, 6531.0, 6658.0, 6787.0, 6918.0, 7052.0, 7188.0, 7326.0, 7467.0, 7610.0, 7756.0, 7904.0, 8055.0, 8208.0, 8364.0, 8523.0, 8685.0, 8850.0, 9019.0, 9191.0, 9366.0, 9544.0, 9726.0, 9911.0, 10100.0, 10292.0, 10488.0, 10688.0, 10892.0, 11100.0, 11312.0, 11528.0, 11748.0, 11972.0, 12200.0, 12432.0, 12669.0, 12910.0, 13156.0",
                           "m_Kinematic_yStar_Binning"               : "-0.7, -0.63, -0.56, -0.49, -0.42, -0.35, -0.28, -0.21, -0.14, -0.07, 0.0, 0.07, 0.14, 0.21, 0.28, 0.35, 0.42, 0.49, 0.56, 0.63. 0.7",
                           "m_Kinematic_deltaPhi_Binning"            : "0.001, 0.14, 0.28, 0.42, 0.56, 0.7, 0. 84, 0.98, 1.12, 1.26, 1.4, 1.54, 1.68, 1.82, 1.96, 2.1, 2.24, 2.38, 2.52, 2.66, 2.8, 2.94, 3.08, 3.22, 3.36, 3.5 ",
                           "m_Kinematic_pTBalance_Binning"           : "0.001, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0",
                           "m_Kinematic_HT_Binning"                  : "15. ,20. ,25. ,35. ,45. ,55. ,70. ,85. ,100. ,116. ,134. ,152. ,172. ,194. ,216. ,240. ,264. ,290. ,318. ,346.,376.,408.,442.,478.,516.,556.,598.,642.,688.,736.,786.,838.,894.,952.,1012.,1076.,1162.,1250.,1310.,1420.,1530.,1750.,1992.,2242.,2500.,2850.,3200.,3600.,4000.,4600.",
                           "m_Kinematic_MHT_Binning"                 : "15. ,20. ,25. ,35. ,45. ,55. ,70. ,85. ,100. ,116. ,134. ,152. ,172. ,194. ,216. ,240. ,264. ,290. ,318. ,346.,376.,408.,442.,478.,516.,556.,598.,642.,688.,736.,786.,838.,894.,952.,1012.,1076.,1162.,1250.,1310.,1420.,1530.,1750.,1992.,2242.,2500.,2850.,3200.,3600.,4000.,4600.",
                           "m_Kinematic_MHTPhi_Binning"              : "-3.5, -3.36, -3.22, -3.08, -2.94, -2.8, -2.66, -2.52, -2.38, -2.24, -2.1, -1.96, -1.82, -1.68, -1.54, -1.4, -1.26, -1.12, -0.98, -0.84, -0.7, -0.56, -0.42, -0.28, -0.14, 0.0, 0.14, 0.28, 0.42, 0.56, 0.7, 0. 84, 0.98, 1.12, 1.26, 1.4, 1.54, 1.68, 1.82, 1.96, 2.1, 2.24, 2.38, 2.52, 2.66, 2.8, 2.94, 3.08, 3.22, 3.36, 3.5",
                           "m_Kinematic_EMFrac_Binning"              : "0.001, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0",
                           "m_Kinematic_HECFrac_Binning"             : "0.001, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0",
                           "m_Kinematic_FracSamplingMax_Binning"     : "0.001, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0",
#
# ============= 2. Response ==========================================================================================================================
#
# Switch Matching on or off
#
                           "m_doMatching"                            : True,
#
# Set the maximum eta-phi-distance allowed for matching (DeltaRMax):
#
                           "m_DeltaRMax"                             : 0.4,
#
# Enable the generation of Pt Response Plots for offline, trigger and/or truth jets
#
                           "m_doTriggerOfflineResponse"              : False,
                           "m_doOfflineTruthResponse"                : False,
                           "m_doTriggerTruthResponse"                : False,
#
# Enable the generation of Mjj Response Plots for offline, trigger and/or truth jets
#
                           "m_doMjjResponseTrigVsOff"                : False,
                           "m_doMjjResponseOffVsTruth"               : False,
                           "m_doMjjResponseTrigVsTruth"              : False,
#
# Binning of Response Plots themselves:
#
                           "m_Response_BinNumbers"                   : 100,
                           "m_Response_MinBin"                       : 0.5,
                           "m_Response_MaxBin"                       : 1.5,
#
# Pt/Mjj and Eta Binning. Here you define the binning in pt and eta creating a histogram for each bin in the pt-eta-space:
# elements have to be separated by kommas ","
#
                           "m_PtResponse_ptBinning"                  : "30.0, 50.0, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510",
                           "m_PtResponse_etaBinning"                 : "-4.9,-3.1,-2.8,-2.1,-1.8,-1.5,-1.2,-0.8,0.0,0.8,1.2,1.5,1.8,2.1,2.8,3.1,4.9",
#
                           "m_MjjResponse_mjjBinning"                : "50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 203.0, 216.0, 229.0, 243.0, 257.0, 272.0, 287.0, 303.0, 319.0, 335.0, 352.0, 369.0, 387.0, 405.0, 424.0, 443.0, 462.0, 482.0, 502.0, 523.0, 544.0, 566.0, 588.0, 611.0, 634.0, 657.0, 681.0, 705.0, 730.0, 755.0, 781.0, 807.0, 834.0, 861.0, 889.0, 917.0, 946.0, 976.0, 1006.0, 1037.0, 1068.0, 1100.0, 1133.0, 1166.0, 1200.0, 1234.0, 1269.0, 1305.0, 1341.0, 1378.0, 1416.0, 1454.0, 1493.0, 1533.0, 1573.0, 1614.0, 1656.0, 1698.0, 1741.0, 1785.0, 1830.0, 1875.0, 1921.0, 1968.0, 2016.0, 2065.0, 2114.0, 2164.0, 2215.0, 2267.0, 2320.0, 2374.0, 2429.0, 2485.0, 2542.0, 2600.0, 2659.0, 2719.0, 2780.0, 2842.0, 2905.0, 2969.0, 3034.0, 3100.0, 3167.0, 3235.0, 3305.0, 3376.0, 3448.0, 3521.0, 3596.0, 3672.0, 3749.0, 3827.0, 3907.0, 3988.0, 4070.0, 4154.0, 4239.0, 4326.0, 4414.0, 4504.0, 4595.0, 4688.0, 4782.0, 4878.0, 4975.0, 5074.0, 5175.0, 5277.0, 5381.0, 5487.0, 5595.0, 5705.0, 5817.0, 5931.0, 6047.0, 6165.0, 6285.0, 6407.0, 6531.0, 6658.0, 6787.0, 6918.0, 7052.0, 7188.0, 7326.0, 7467.0, 7610.0, 7756.0, 7904.0, 8055.0, 8208.0, 8364.0, 8523.0, 8685.0, 8850.0, 9019.0, 9191.0, 9366.0, 9544.0, 9726.0, 9911.0, 10100.0, 10292.0, 10488.0, 10688.0, 10892.0, 11100.0, 11312.0, 11528.0, 11748.0, 11972.0, 12200.0, 12432.0, 12669.0, 12910.0, 13156.0",
                           "m_MjjResponse_etaBinning"                : "-4.9, 4.9",
#
# ============= 3. TurnOns (Trigger Efficiencies) ===================================================================================================
#
#
# Enable the generation of turnons
#
                           "m_doTurnOns"                             : True, 
#
# Select the STRATEGY of turnOn generation.
# The available options are the Trigger Decision Tool (TDT), Emulation or Trigger Before Preascale. Note that more than one strategy can be selected
# and will be created independently of each other
#
                           "m_useTriggerDecisionTool"                : True,
                           "m_useEmulation"                          : False,
                           "m_useTriggerBeforePraescale"             : False,
#
# Select triggers which turnon plots you would like to generate.
# The syntax for turnon plots is 'probe-trigger'-'reference-trigger', e.g. HLT_j25-HLT_j15
# In brackets, the strategy of turnon generation can be selected specifically for each trigger, e.g. (tep)
# The options for the brackets are: 't'= TDT, 'e'= Emulation, 'p' = TBP. Without any bracket, all three options are selected by default.
#
                           "m_TurnOnName"                            : "HLT_j25/HLT_j15 (tp);  HLT_j60/HLT_j25; HLT_j110/HLT_j85; HLT_j360/HLT_j260;  HLT_j380/HLT_j260; HLT_j400/  HLT_j260",
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
#
# Binning of pt, HT and TLA turnons:
#
                           "m_PtTurnon_BinNumbers"                   : 1200,
                           "m_PtTurnon_MinBin"                       : 0.0,
                           "m_PtTurnon_MaxBin"                       : 1200.0,
#
                           "m_HtTurnon_BinNumbers"                   : 1200,
                           "m_HtTurnon_MinBin"                       : 0.0,
                           "m_HtTurnon_MaxBin"                       : 1200.0,
#
                           "m_MjjTurnon_BinNumbers"                  : 1200,
                           "m_MjjTurnon_MinBin"                      : 0.0,
                           "m_MjjTurnon_MaxBin"                      : 1200.0,
#
# CLEANING                           "m_doCleaning"                            : False,
                           } )

