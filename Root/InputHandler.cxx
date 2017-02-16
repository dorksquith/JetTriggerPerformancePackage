/*--------------------------------------------------------------------
InputHandler.cxx

The InputHandler is JTPP's main function where all wanted histogram classes
are declared and their methods are used.

created by Edgar Kellermann (edgar.kellermann@cern.ch)
--------------------------------------------------------------------*/
#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include "EventLoop/OutputStream.h"
#include "AthContainers/ConstDataVector.h"

#include "xAODTracking/VertexContainer.h"
#include <xAODJet/JetContainer.h>
#include "xAODEventInfo/EventInfo.h"
#include <JetTriggerPerformancePackage/InputHandler.h>
#include <xAODAnaHelpers/HelperFunctions.h>
#include <xAODAnaHelpers/tools/ReturnCheck.h>

#include "TFile.h"
#include "TKey.h"
#include "TLorentzVector.h"
#include "TEnv.h"
#include "TSystem.h"

#include <utility>
#include <iostream>
#include <fstream>


using namespace std;

// this is needed to distribute the algorithm to the workers
ClassImp(InputHandler)

InputHandler :: InputHandler () :
    //constructor for variables
    m_applyGRL(true),
    m_GRLxml("$ROOTCOREBIN/data/ZprimeDM/data15_13TeV.periodAllYear_DetStatus-v64-pro19_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml"),
    m_branch_offlineJet("jet_X"),
    m_branch_triggerJet("trigJet_X"),
    m_branch_truthJet("truthJet_X"),
    m_branch_LVL1JetROI("LVL1JetROIs_X"),
    m_TriggerName("STUDYALL"),
    m_debug(false),
    m_debugInExecute(false),
    m_coutPassedTriggers(false),
    m_isData(true),
    m_checkLArError(false),
    m_doOfflineJetKinematics(false),
    m_doTriggerJetKinematics(false),
    m_doTruthJetKinematics(false),
    m_doLeadSubleadThirdJetKinematics(false),
    m_doNthJetKinematics(false),
    m_doKinematicsInBinsOfEta(false),
    m_doEtainBinsOfPt(false),
    m_doMjj(false),
    m_doM23(false),
    m_doyStar(false),
    m_doDeltaPhi(false),
    m_doPTBalance(false),
    m_doHT(false),
    m_doMHT(false),
    m_doMHTPhi(false),
    m_doEMFrac(false),
    m_doHECFrac(false),
    m_doFracSamplingMax(false),
    m_calculateMjj(false),
    m_calculateM23(false),
    m_calculateYStar(false),
    m_Kinematic_E_Binning("15. ,20. ,25. ,35. ,45. ,55. ,70. ,85. ,100. ,116. ,134. ,152. ,172. ,194. ,216. ,240. ,264. ,290. ,318. ,346.,376.,408.,442.,478.,516.,556.,598.,642.,688.,736.,786.,838.,894.,952.,1012.,1076.,1162.,1250.,1310.,1420.,1530.,1750.,1992.,2242.,2500.,2850.,3200.,3600.,4000.,4600."),
    m_Kinematic_pt_Binning("15. ,20. ,25. ,35. ,45. ,55. ,70. ,85. ,100. ,116. ,134. ,152. ,172. ,194. ,216. ,240. ,264. ,290. ,318. ,346.,376.,408.,442.,478.,516.,556.,598.,642.,688.,736.,786.,838.,894.,952.,1012.,1076.,1162.,1250.,1310.,1420.,1530.,1750.,1992.,2242.,2500.,2850.,3200.,3600.,4000.,4600."),
    m_Kinematic_phi_Binning("-3.5, -3.36, -3.22, -3.08, -2.94, -2.8, -2.66, -2.52, -2.38, -2.24, -2.1, -1.96, -1.82, -1.68, -1.54, -1.4, -1.26, -1.12, -0.98, -0.84, -0.7, -0.56, -0.42, -0.28, -0.14, 0.0, 0.14, 0.28, 0.42, 0.56, 0.7, 0. 84, 0.98, 1.12, 1.26, 1.4, 1.54, 1.68, 1.82, 1.96, 2.1, 2.24, 2.38, 2.52, 2.66, 2.8, 2.94, 3.08, 3.22, 3.36, 3.5"),
    m_Kinematic_eta_Binning("-4.9,-3.1,-2.8,-2.1,-1.8,-1.5,-1.2,-0.8,0.0,0.8,1.2,1.5,1.8,2.1,2.8,3.1,4.9"),
    m_Kinematic_mjj_Binning("50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 203.0, 216.0, 229.0, 243.0, 257.0, 272.0, 287.0, 303.0, 319.0, 335.0, 352.0, 369.0, 387.0, 405.0, 424.0, 443.0, 462.0, 482.0, 502.0, 523.0, 544.0, 566.0, 588.0, 611.0, 634.0, 657.0, 681.0, 705.0, 730.0, 755.0, 781.0, 807.0, 834.0, 861.0, 889.0, 917.0, 946.0, 976.0, 1006.0, 1037.0, 1068.0, 1100.0, 1133.0, 1166.0, 1200.0, 1234.0, 1269.0, 1305.0, 1341.0, 1378.0, 1416.0, 1454.0, 1493.0, 1533.0, 1573.0, 1614.0, 1656.0, 1698.0, 1741.0, 1785.0, 1830.0, 1875.0, 1921.0, 1968.0, 2016.0, 2065.0, 2114.0, 2164.0, 2215.0, 2267.0, 2320.0, 2374.0, 2429.0, 2485.0, 2542.0, 2600.0, 2659.0, 2719.0, 2780.0, 2842.0, 2905.0, 2969.0, 3034.0, 3100.0, 3167.0, 3235.0, 3305.0, 3376.0, 3448.0, 3521.0, 3596.0, 3672.0, 3749.0, 3827.0, 3907.0, 3988.0, 4070.0, 4154.0, 4239.0, 4326.0, 4414.0, 4504.0, 4595.0, 4688.0, 4782.0, 4878.0, 4975.0, 5074.0, 5175.0, 5277.0, 5381.0, 5487.0, 5595.0, 5705.0, 5817.0, 5931.0, 6047.0, 6165.0, 6285.0, 6407.0, 6531.0, 6658.0, 6787.0, 6918.0, 7052.0, 7188.0, 7326.0, 7467.0, 7610.0, 7756.0, 7904.0, 8055.0, 8208.0, 8364.0, 8523.0, 8685.0, 8850.0, 9019.0, 9191.0, 9366.0, 9544.0, 9726.0, 9911.0, 10100.0, 10292.0, 10488.0, 10688.0, 10892.0, 11100.0, 11312.0, 11528.0, 11748.0, 11972.0, 12200.0, 12432.0, 12669.0, 12910.0, 13156.0"),
    m_Kinematic_m23_Binning("50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 203.0, 216.0, 229.0, 243.0, 257.0, 272.0, 287.0, 303.0, 319.0, 335.0, 352.0, 369.0, 387.0, 405.0, 424.0, 443.0, 462.0, 482.0, 502.0, 523.0, 544.0, 566.0, 588.0, 611.0, 634.0, 657.0, 681.0, 705.0, 730.0, 755.0, 781.0, 807.0, 834.0, 861.0, 889.0, 917.0, 946.0, 976.0, 1006.0, 1037.0, 1068.0, 1100.0, 1133.0, 1166.0, 1200.0, 1234.0, 1269.0, 1305.0, 1341.0, 1378.0, 1416.0, 1454.0, 1493.0, 1533.0, 1573.0, 1614.0, 1656.0, 1698.0, 1741.0, 1785.0, 1830.0, 1875.0, 1921.0, 1968.0, 2016.0, 2065.0, 2114.0, 2164.0, 2215.0, 2267.0, 2320.0, 2374.0, 2429.0, 2485.0, 2542.0, 2600.0, 2659.0, 2719.0, 2780.0, 2842.0, 2905.0, 2969.0, 3034.0, 3100.0, 3167.0, 3235.0, 3305.0, 3376.0, 3448.0, 3521.0, 3596.0, 3672.0, 3749.0, 3827.0, 3907.0, 3988.0, 4070.0, 4154.0, 4239.0, 4326.0, 4414.0, 4504.0, 4595.0, 4688.0, 4782.0, 4878.0, 4975.0, 5074.0, 5175.0, 5277.0, 5381.0, 5487.0, 5595.0, 5705.0, 5817.0, 5931.0, 6047.0, 6165.0, 6285.0, 6407.0, 6531.0, 6658.0, 6787.0, 6918.0, 7052.0, 7188.0, 7326.0, 7467.0, 7610.0, 7756.0, 7904.0, 8055.0, 8208.0, 8364.0, 8523.0, 8685.0, 8850.0, 9019.0, 9191.0, 9366.0, 9544.0, 9726.0, 9911.0, 10100.0, 10292.0, 10488.0, 10688.0, 10892.0, 11100.0, 11312.0, 11528.0, 11748.0, 11972.0, 12200.0, 12432.0, 12669.0, 12910.0, 13156.0"),
    m_Kinematic_yStar_Binning("-0.7, -0.63, -0.56, -0.49, -0.42, -0.35, -0.28, -0.21, -0.14, -0.07, 0.0, 0.07, 0.14, 0.21, 0.28, 0.35, 0.42, 0.49, 0.56, 0.63. 0.7"),
    m_Kinematic_deltaPhi_Binning("0.001, 0.14, 0.28, 0.42, 0.56, 0.7, 0. 84, 0.98, 1.12, 1.26, 1.4, 1.54, 1.68, 1.82, 1.96, 2.1, 2.24, 2.38, 2.52, 2.66, 2.8, 2.94, 3.08, 3.22, 3.36, 3.5 "),
    m_Kinematic_pTBalance_Binning("0.001, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0"),
    m_Kinematic_HT_Binning("15. ,20. ,25. ,35. ,45. ,55. ,70. ,85. ,100. ,116. ,134. ,152. ,172. ,194. ,216. ,240. ,264. ,290. ,318. ,346.,376.,408.,442.,478.,516.,556.,598.,642.,688.,736.,786.,838.,894.,952.,1012.,1076.,1162.,1250.,1310.,1420.,1530.,1750.,1992.,2242.,2500.,2850.,3200.,3600.,4000.,4600."),
    m_Kinematic_MHT_Binning("15. ,20. ,25. ,35. ,45. ,55. ,70. ,85. ,100. ,116. ,134. ,152. ,172. ,194. ,216. ,240. ,264. ,290. ,318. ,346.,376.,408.,442.,478.,516.,556.,598.,642.,688.,736.,786.,838.,894.,952.,1012.,1076.,1162.,1250.,1310.,1420.,1530.,1750.,1992.,2242.,2500.,2850.,3200.,3600.,4000.,4600."),
    m_Kinematic_MHTPhi_Binning("-3.5, -3.36, -3.22, -3.08, -2.94, -2.8, -2.66, -2.52, -2.38, -2.24, -2.1, -1.96, -1.82, -1.68, -1.54, -1.4, -1.26, -1.12, -0.98, -0.84, -0.7, -0.56, -0.42, -0.28, -0.14, 0.0, 0.14, 0.28, 0.42, 0.56, 0.7, 0. 84, 0.98, 1.12, 1.26, 1.4, 1.54, 1.68, 1.82, 1.96, 2.1, 2.24, 2.38, 2.52, 2.66, 2.8, 2.94, 3.08, 3.22, 3.36, 3.5"),
    m_Kinematic_EMFrac_Binning("0.001, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0"),
    m_Kinematic_HECFrac_Binning("0.001, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0"),
    m_Kinematic_FracSamplingMax_Binning("0.001, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0"),
    m_doMatching(false),
    m_DeltaRMax(0.4),
    m_doOfflineTruthResponse(false),
    m_doTriggerTruthResponse(false),
    m_doTriggerOfflineResponse(false),
    m_Response_BinNumbers(100),
    m_Response_MinBin(0.5),
    m_Response_MaxBin(1.5),
    m_doMjjResponseTrigVsOff(false),
    m_doMjjResponseOffVsTruth(false),
    m_doMjjResponseTrigVsTruth(false),
    m_PtResponse_ptBinning("30.0, 50.0, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510"),
    m_PtResponse_etaBinning("-4.9,-3.1,-2.8,-2.1,-1.8,-1.5,-1.2,-0.8,0.0,0.8,1.2,1.5,1.8,2.1,2.8,3.1,4.9"),
    m_MjjResponse_mjjBinning("50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 203.0, 216.0, 229.0, 243.0, 257.0, 272.0, 287.0, 303.0, 319.0, 335.0, 352.0, 369.0, 387.0, 405.0, 424.0, 443.0, 462.0, 482.0, 502.0, 523.0, 544.0, 566.0, 588.0, 611.0, 634.0, 657.0, 681.0, 705.0, 730.0, 755.0, 781.0, 807.0, 834.0, 861.0, 889.0, 917.0, 946.0, 976.0, 1006.0, 1037.0, 1068.0, 1100.0, 1133.0, 1166.0, 1200.0, 1234.0, 1269.0, 1305.0, 1341.0, 1378.0, 1416.0, 1454.0, 1493.0, 1533.0, 1573.0, 1614.0, 1656.0, 1698.0, 1741.0, 1785.0, 1830.0, 1875.0, 1921.0, 1968.0, 2016.0, 2065.0, 2114.0, 2164.0, 2215.0, 2267.0, 2320.0, 2374.0, 2429.0, 2485.0, 2542.0, 2600.0, 2659.0, 2719.0, 2780.0, 2842.0, 2905.0, 2969.0, 3034.0, 3100.0, 3167.0, 3235.0, 3305.0, 3376.0, 3448.0, 3521.0, 3596.0, 3672.0, 3749.0, 3827.0, 3907.0, 3988.0, 4070.0, 4154.0, 4239.0, 4326.0, 4414.0, 4504.0, 4595.0, 4688.0, 4782.0, 4878.0, 4975.0, 5074.0, 5175.0, 5277.0, 5381.0, 5487.0, 5595.0, 5705.0, 5817.0, 5931.0, 6047.0, 6165.0, 6285.0, 6407.0, 6531.0, 6658.0, 6787.0, 6918.0, 7052.0, 7188.0, 7326.0, 7467.0, 7610.0, 7756.0, 7904.0, 8055.0, 8208.0, 8364.0, 8523.0, 8685.0, 8850.0, 9019.0, 9191.0, 9366.0, 9544.0, 9726.0, 9911.0, 10100.0, 10292.0, 10488.0, 10688.0, 10892.0, 11100.0, 11312.0, 11528.0, 11748.0, 11972.0, 12200.0, 12432.0, 12669.0, 12910.0, 13156.0"),
    m_MjjResponse_etaBinning("-4.9, 4.9"),
    m_doTurnOns(false),
    m_TurnOnName(""),
    m_useTriggerDecisionTool(false),
    m_useEmulation(false),
    m_useTriggerBeforePraescale(false),
    m_PtTurnon_BinNumbers(1200),
    m_PtTurnon_MinBin(0.0),
    m_PtTurnon_MaxBin(1200.0),
    m_HtTurnon_BinNumbers(1200),
    m_HtTurnon_MinBin(0.0),
    m_HtTurnon_MaxBin(1200.0),
    m_MjjTurnon_BinNumbers(1200),
    m_MjjTurnon_MinBin(0.0),
    m_MjjTurnon_MaxBin(1200.0),
    m_HLT_cen_etaMin(0.0), // DEFAULT ETA REGIONS AND TIMING FOR TURNON EVENT SELECTION!
    m_HLT_cen_etaMax(2.8), //
    m_HLT_fwd_etaMin(3.6), //
    m_HLT_fwd_etaMax(99.0), //
    m_L1_cen_etaMin(0.0),  //
    m_L1_cen_etaMax(2.6),  //
    m_L1_fwd_etaMin(3.6),  //
    m_L1_fwd_etaMax(99.0),  //
    m_TurnonCut_Timing(10.0), //
    m_doCuts(false),
    m_cutStringKinematics(""),
    m_cutStringTurnons(""),
    m_doCleaning(false),
    m_doOnlyThisNumberOfEvents(-1),
    cutH(nullptr),
    triggerhistos(nullptr),
    trigEffic(nullptr),
    TD(nullptr),
    ED_jet(nullptr),
    ED_trigJet(nullptr),
    ED_truthJet(nullptr),
    L1D(nullptr),
    anaHandler(nullptr),
    logwr(nullptr),
    out_trigger(nullptr),
    out_turnOns(nullptr),
    out_logfile(nullptr),
    m_testArray(nullptr)
{
  Info("InputHandler()", "Calling constructor");
}

// Declare and construct all needed classes
EL::StatusCode  InputHandler :: configure ()
{
  if(m_debug) Info("InputHandler()", "Calling InputHandler :: configure");  

  // Declare and construct a ConfigStatus class
  CS = new ConfigStatus(m_debug,
                        m_debugInExecute,
                        m_TriggerName,
                        m_isData,
                        m_checkLArError,
                        m_doOfflineJetKinematics,
                        m_doTriggerJetKinematics,
                        m_doTruthJetKinematics,
                        m_doLeadSubleadThirdJetKinematics,
                        m_doNthJetKinematics,
                        m_doKinematicsInBinsOfEta,
                        m_doEtainBinsOfPt,
                        m_doMjj,
                        m_doM23,
                        m_doyStar,
                        m_doDeltaPhi,
                        m_doPTBalance,
                        m_doHT,
                        m_doMHT,
                        m_doMHTPhi,
                        m_doEMFrac,
                        m_doHECFrac,
                        m_doFracSamplingMax,
                        m_calculateMjj,
                        m_calculateM23,
                        m_calculateYStar,
                        m_Kinematic_E_Binning,
                        m_Kinematic_pt_Binning,
                        m_Kinematic_phi_Binning,
                        m_Kinematic_eta_Binning,
                        m_Kinematic_mjj_Binning,
                        m_Kinematic_m23_Binning,
                        m_Kinematic_yStar_Binning,
                        m_Kinematic_deltaPhi_Binning,
                        m_Kinematic_pTBalance_Binning,
                        m_Kinematic_HT_Binning,
                        m_Kinematic_MHT_Binning,
                        m_Kinematic_MHTPhi_Binning,
                        m_Kinematic_EMFrac_Binning,
                        m_Kinematic_HECFrac_Binning,
                        m_Kinematic_FracSamplingMax_Binning,
                        m_doMatching,
                        m_doOfflineTruthResponse,
                        m_doTriggerTruthResponse,
                        m_doTriggerOfflineResponse,
                        m_doMjjResponseTrigVsOff,
                        m_doMjjResponseOffVsTruth,
                        m_doMjjResponseTrigVsTruth,
                        m_Response_BinNumbers,
                        m_Response_MinBin,
                        m_Response_MaxBin,
                        m_PtResponse_ptBinning,
                        m_PtResponse_etaBinning,
                        m_MjjResponse_mjjBinning,
                        m_MjjResponse_etaBinning,
                        m_doTurnOns,
                        m_TurnOnName,
                        m_useTriggerDecisionTool,
                        m_useEmulation,
                        m_useTriggerBeforePraescale,
                        m_PtTurnon_BinNumbers,
                        m_PtTurnon_MinBin,
                        m_PtTurnon_MaxBin,
                        m_HtTurnon_BinNumbers,
                        m_HtTurnon_MinBin,
                        m_HtTurnon_MaxBin,
                        m_MjjTurnon_BinNumbers,
                        m_MjjTurnon_MinBin,
                        m_MjjTurnon_MaxBin,
                        m_HLT_cen_etaMin,
                        m_HLT_cen_etaMax,
                        m_HLT_fwd_etaMin,
                        m_HLT_fwd_etaMax,
                        m_L1_cen_etaMin,
                        m_L1_cen_etaMax,
                        m_L1_fwd_etaMin,
                        m_L1_fwd_etaMax,
                        m_TurnonCut_Timing,
                        m_doCuts,
                        m_cutStringKinematics,
                        m_cutStringTurnons,
                        m_doCleaning,
                        m_doOnlyThisNumberOfEvents);

  // Declare and construct CutHandler for kinematic plots
  if (CS->doCuts){
    std::cout << "\n==== Kinematics Cut Selection ====" << std::endl;
    cutH = new CutHandler(CS->cutStringKinematics);
  }
  else cutH = new CutHandler();

  // Declare and construct a TriggerData class
  TD = new TriggerData(m_TriggerName, m_TurnOnName, CS);

  // Logfiles for JTPPPlotter.py
  out_trigger = new ofstream("logfile_trigger.txt");
  out_turnOns = new ofstream("logfile_turnOns.txt");
  out_logfile = new ofstream("logfileJTPP.txt");

  // Declare and construct LogWriter
  logwr = new LogWriter("logJTPP.txt");

  // eventCounter variable in order to keep track of the number if events
  eventCounter = 0;

  // Declare and construct ToolsJTPP();
  ToolsJTPP* myTools = new ToolsJTPP();

  // Cut branch names and fill variables
  m_branch_offlineJet_front = myTools->splitString(m_branch_offlineJet,"X",0);
  m_branch_offlineJet_back  = myTools->splitString(m_branch_offlineJet,"X",1);

  m_branch_triggerJet_front = myTools->splitString(m_branch_triggerJet,"X",0);
  m_branch_triggerJet_back  = myTools->splitString(m_branch_triggerJet,"X",1);

  m_branch_truthJet_front = myTools->splitString(m_branch_truthJet,"X",0);
  m_branch_truthJet_back  = myTools->splitString(m_branch_truthJet,"X",1);

  m_branch_LVL1JetROI_front = myTools->splitString(m_branch_LVL1JetROI,"X",0);
  m_branch_LVL1JetROI_back  = myTools->splitString(m_branch_LVL1JetROI,"X",1);

  // write all selected triggers into logfiles
  *out_logfile << "\n===Triggers===" << endl;
  for (unsigned int i=0; i<TD->config_triggerName.size(); i++){
      *out_trigger << TD->config_triggerName.at(i) << endl;
      *out_logfile << TD->config_triggerName.at(i) << endl;
  }

  // write all selected triggers into logfiles
  for (unsigned int i=0; i<TD->probe_triggerName.size(); i++){
      *out_turnOns << TD->probe_triggerName.at(i) << "-" << TD->ref_triggerName.at(i) << endl;
  }

  // Declare a TriggerEfficiencyMatrix class (containing all Turn-ons) and book all wanted histograms
  if (CS->doTurnOns){
      trigEffic = new TriggerEfficiencyMatrix("effic", "TurnOns", CS);
      trigEffic->BookAll(TD, CS, wk());
  }

  // Declare all TriggerHistoPack classes (one for each trigger) and book all histograms
  for (unsigned int i=0; i<TD->config_triggerName.size(); i++){

      TriggerHistoPack* triggerhisto = new TriggerHistoPack(TD->config_triggerName.at(i), TD->config_triggerName.at(i), TD->config_nthJet.at(i), CS);
      triggerhisto->BookAll(wk());
      m_map_triggerhistos[TD->config_triggerName.at(i)] = triggerhisto;
  }

  // Declare and construct EventData classes for all jet types, i.e. offline jets, trigger jets, truth jets, matched truth jets and L1 jets
  ED_jet = new EventData("jet", CS);
  ED_trigJet = new EventData("trigJet", CS);
  ED_truthJet = new EventData("truthJet", CS);
  L1D = new L1Data("L1Jet");

  // Declare the AnalysisHandler class
  anaHandler = new AnalysisHandler();

  // Fill logfile with more informations needed in JTPPPlotter
  // 1. Leading, subleading, third jet
  *out_logfile << "===LeadSubleadThird===" << endl;
  *out_logfile << "" << endl; // empty space for all jets
  if (m_doLeadSubleadThirdJetKinematics){
      *out_logfile << "lead" << endl;
      *out_logfile << "sublead" << endl;
      *out_logfile << "third" << endl;
  }

  // Fill ED->isGood vector at the first time:
  if(m_debug) std::cout << "Initialising ED->isGood..." << std::endl;
  for (unsigned int i=0; i< 999; i++){
      ED_jet->isGood.push_back(true);
      ED_trigJet->isGood.push_back(true);
      ED_truthJet->isGood.push_back(true);
  }

  if(m_debug) Info("InputHandler()", "Finishing InputHandler :: configure");

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode InputHandler :: setupJob (EL::Job& job)
{
  if(m_debug) std::cout << "EDDI: Calling InputHandler :: setupJob" << std::endl;
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode InputHandler :: histInitialize ()
{
  if(m_debug) std::cout << "EDDI: Calling InputHandler :: histInitialize" << std::endl;
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  Info("histInitialize()", "Calling histInitialize \n");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode InputHandler :: fileExecute ()
{
  if(m_debug) std::cout << "EDDI: Calling InputHandler :: fileExecute" << std::endl;
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}


// Here we open the file(s) and set the needed branches
EL::StatusCode InputHandler :: changeInput (bool firstFile)
{
  if(m_debug) Info("InputHandler()", "Calling InputHandler :: changeInput");


  // Opening file
  if(firstFile){
    if ( this->configure() == EL::StatusCode::FAILURE ) {
      Error("initialize()", "Failed to properly configure. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    TFile* inputFile = wk()->inputFile();
    TIter next(inputFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
      std::string keyName = key->GetName();

      std::size_t found = keyName.find("cutflow");
      bool foundCutFlow = (found!=std::string::npos);

      found = keyName.find("weighted");
      bool foundWeighted = (found!=std::string::npos);

    }//over Keys
    if(m_debug) std::cout << "looking at cutflow..." << std::endl;

    //for MC events: open cutflow histo and obtain numberOfEntries
    if(!CS->isData){
        numberOfEntries = 1; // For some nutples, cutflow does not exist
    }

  }// first file

  TTree *tree = wk()->tree();
  if(m_debug) Info("InputHandler()", "Setting default branches");

  tree->SetBranchStatus ("*", 0);

  tree->SetBranchStatus  ("passedTriggers",    1);
  tree->SetBranchAddress ("passedTriggers",    &(TD->event_passedTriggers));

  if ((CS->doTurnOns) && CS->useTriggerBeforePraescale){

    tree->SetBranchStatus  ("triggerPrescales",    1);
    tree->SetBranchAddress ("triggerPrescales",  &(TD->event_triggerPrescales));

    tree->SetBranchStatus  ("isPassBits",    1);
    tree->SetBranchAddress ("isPassBits",  &(TD->event_isPassBits));

    tree->SetBranchStatus  ("isPassBitsNames",    1);
    tree->SetBranchAddress ("isPassBitsNames",  &(TD->event_isPassBitsNames));
  }

  tree->SetBranchStatus  ("runNumber",    1);
  tree->SetBranchAddress ("runNumber",    &m_runNumber);

  tree->SetBranchStatus  ("eventNumber",    1);
  tree->SetBranchAddress ("eventNumber",    &m_eventNumber);

  tree->SetBranchStatus  ("lumiBlock",    1);
  tree->SetBranchAddress ("lumiBlock",    &m_lumiBlock);

  tree->SetBranchStatus  ("weight",    1);
  tree->SetBranchAddress ("weight",    &m_weight);


  // Data specific observables
  if ((CS->isData) && (CS->checkLArError)){
      tree->SetBranchStatus  ("LArError", 1);
      tree->SetBranchAddress ("LArError", &m_LArError);
  }

  // Offline jets
  if ((CS->doOfflineJetKinematics)||(CS->doOfflineTruthResponse)||(CS->doTriggerOfflineResponse)|| (CS->doMjjResponseOffVsTruth) || (CS->doMjjResponseTrigVsOff) || (CS->doTurnOns)){

      if(m_debug) Info("InputHandler()", "Setting offline branches");

      tree->SetBranchStatus  ((m_branch_offlineJet_front + "pt" + m_branch_offlineJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet_front + "pt" + m_branch_offlineJet_back).c_str(), &(ED_jet->pt));

      tree->SetBranchStatus  ((m_branch_offlineJet_front + "eta" + m_branch_offlineJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet_front + "eta" + m_branch_offlineJet_back).c_str(), &(ED_jet->eta));

      tree->SetBranchStatus  ((m_branch_offlineJet_front + "phi" + m_branch_offlineJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet_front + "phi" + m_branch_offlineJet_back).c_str(), &(ED_jet->phi));

      tree->SetBranchStatus  ((m_branch_offlineJet_front + "E" + m_branch_offlineJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet_front + "E" + m_branch_offlineJet_back).c_str(), &(ED_jet->E));

      tree->SetBranchStatus  ((m_branch_offlineJet_front + "mjj" + m_branch_offlineJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet_front + "mjj" + m_branch_offlineJet_back).c_str(), &(ED_jet->mjj));

      tree->SetBranchStatus  ((m_branch_offlineJet_front + "m23" + m_branch_offlineJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet_front + "m23" + m_branch_offlineJet_back).c_str(), &(ED_jet->m23));

      if ((CS->isData) && (CS->doCleaning)){
	  tree->SetBranchStatus  ((m_branch_offlineJet_front + "clean_passLooseBad" + m_branch_offlineJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_offlineJet_front + "clean_passLooseBad" + m_branch_offlineJet_back).c_str(), &(ED_jet->passedCleaning));
      }

      if(CS->doyStar){
	  tree->SetBranchStatus  ((m_branch_offlineJet_front + "yStar" + m_branch_offlineJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_offlineJet_front + "yStar" + m_branch_offlineJet_back).c_str(), &(ED_jet->yStar));
      }
      if(CS->doDeltaPhi){
	  tree->SetBranchStatus  ((m_branch_offlineJet_front + "deltaPhi" + m_branch_offlineJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_offlineJet_front + "deltaPhi" + m_branch_offlineJet_back).c_str(), &(ED_jet->deltaPhi));
      }
      if(CS->doPTBalance){
	  tree->SetBranchStatus  ((m_branch_offlineJet_front + "pTBalance" + m_branch_offlineJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_offlineJet_front + "pTBalance" + m_branch_offlineJet_back).c_str(), &(ED_jet->pTBalance));
      }
      if(CS->doMHT){
	  tree->SetBranchStatus  ((m_branch_offlineJet_front + "MHT" + m_branch_offlineJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_offlineJet_front + "MHT" + m_branch_offlineJet_back).c_str(), &(ED_jet->MHT));
      }
      if(CS->doMHTPhi){
	  tree->SetBranchStatus  ((m_branch_offlineJet_front + "MHTPhi" + m_branch_offlineJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_offlineJet_front + "MHTPhi" + m_branch_offlineJet_back).c_str(), &(ED_jet->MHTPhi));
      }
      if(CS->doEMFrac){
	  tree->SetBranchStatus  ((m_branch_offlineJet_front + "EMFrac" + m_branch_offlineJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_offlineJet_front + "EMFrac" + m_branch_offlineJet_back).c_str(), &(ED_jet->EMFrac));
      }
      if(CS->doHECFrac){
	  tree->SetBranchStatus  ((m_branch_offlineJet_front + "HECFrac" + m_branch_offlineJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_offlineJet_front + "HECFrac" + m_branch_offlineJet_back).c_str(), &(ED_jet->HECFrac));
      }
      if(CS->doFracSamplingMax){
	  tree->SetBranchStatus  ((m_branch_offlineJet_front + "FracSamplingMax" + m_branch_offlineJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_offlineJet_front + "FracSamplingMax" + m_branch_offlineJet_back).c_str(), &(ED_jet->FracSamplingMax));
      }
      if (CS->doTurnOns){
	  tree->SetBranchStatus  ((m_branch_offlineJet_front + "Timing" + m_branch_offlineJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_offlineJet_front + "Timing" + m_branch_offlineJet_back).c_str(), &(ED_jet->Timing));
      }
  }

  // Trigg jets
  if ((CS->doTriggerJetKinematics)||(CS->doTriggerTruthResponse)||(CS->doTriggerOfflineResponse)|| (CS->doMjjResponseTrigVsOff) || (CS->doMjjResponseTrigVsTruth) || ((CS->doTurnOns)&&(CS->useEmulation))){

      if(m_debug) Info("InputHandler()", "Setting trigger branches");

      tree->SetBranchStatus  ((m_branch_triggerJet_front + "pt" + m_branch_triggerJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet_front + "pt" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->pt));

      tree->SetBranchStatus  ((m_branch_triggerJet_front + "eta" + m_branch_triggerJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet_front + "eta" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->eta));

      tree->SetBranchStatus  ((m_branch_triggerJet_front + "phi" + m_branch_triggerJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet_front + "phi" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->phi));

      tree->SetBranchStatus  ((m_branch_triggerJet_front + "E" + m_branch_triggerJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet_front + "E" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->E));

      tree->SetBranchStatus  ((m_branch_triggerJet_front + "mjj" + m_branch_triggerJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet_front + "mjj" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->mjj));

      tree->SetBranchStatus  ((m_branch_triggerJet_front + "m23" + m_branch_triggerJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet_front + "m23" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->m23));

      if ((CS->isData) && (CS->doCleaning)){
	  tree->SetBranchStatus  ((m_branch_triggerJet_front + "_clean_passLooseBad" + m_branch_triggerJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_triggerJet_front + "_clean_passLooseBad" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->passedCleaning));
      }

      if(CS->doyStar){
	  tree->SetBranchStatus  ((m_branch_triggerJet_front + "yStar" + m_branch_triggerJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_triggerJet_front + "yStar" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->yStar));
      }
      if(CS->doDeltaPhi){
	  tree->SetBranchStatus  ((m_branch_triggerJet_front + "deltaPhi" + m_branch_triggerJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_triggerJet_front + "deltaPhi" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->deltaPhi));
      }
      if(CS->doPTBalance){
	  tree->SetBranchStatus  ((m_branch_triggerJet_front + "pTBalance" + m_branch_triggerJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_triggerJet_front + "pTBalance" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->pTBalance));
      }
      if(CS->doMHT){
	  tree->SetBranchStatus  ((m_branch_triggerJet_front + "MHT" + m_branch_triggerJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_triggerJet_front + "MHT" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->MHT));
      }
      if(CS->doMHTPhi){
	  tree->SetBranchStatus  ((m_branch_triggerJet_front + "MHTPhi" + m_branch_triggerJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_triggerJet_front + "MHTPhi" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->MHTPhi));
      }
      if(CS->doEMFrac){
	  tree->SetBranchStatus  ((m_branch_triggerJet_front + "EMFrac" + m_branch_triggerJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_triggerJet_front + "EMFrac" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->EMFrac));
      }
      if(CS->doHECFrac){
	  tree->SetBranchStatus  ((m_branch_triggerJet_front + "HECFrac" + m_branch_triggerJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_triggerJet_front + "HECFrac" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->HECFrac));
      }
      if(CS->doFracSamplingMax){
	  tree->SetBranchStatus  ((m_branch_triggerJet_front + "FracSamplingMax" + m_branch_triggerJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_triggerJet_front + "FracSamplingMax" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->FracSamplingMax));
      }
      if (CS->doTurnOns){
	  tree->SetBranchStatus  ((m_branch_triggerJet_front + "Timing" + m_branch_triggerJet_back).c_str(), 1);
	  tree->SetBranchAddress ((m_branch_triggerJet_front + "Timing" + m_branch_triggerJet_back).c_str(), &(ED_trigJet->Timing));
      }
  }

  // Truth jets
  if ((!CS->isData) &&
      ((CS->doTruthJetKinematics)||(CS->doOfflineTruthResponse)||(CS->doTriggerTruthResponse)|| (CS->doMjjResponseOffVsTruth) || (CS->doMjjResponseTrigVsTruth))){

      if(m_debug) Info("InputHandler()", "Setting truth branches");

      tree->SetBranchStatus  ((m_branch_truthJet_front + "pt" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "pt" + m_branch_truthJet_back).c_str(), &(ED_truthJet->pt));

      tree->SetBranchStatus  ((m_branch_truthJet_front + "eta" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "eta" + m_branch_truthJet_back).c_str(), &(ED_truthJet->eta));

      tree->SetBranchStatus  ((m_branch_truthJet_front + "phi" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "phi" + m_branch_truthJet_back).c_str(), &(ED_truthJet->phi));

      tree->SetBranchStatus  ((m_branch_truthJet_front + "E" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "E" + m_branch_truthJet_back).c_str(), &(ED_truthJet->E));

      tree->SetBranchStatus  ((m_branch_truthJet_front + "mjj" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "mjj" + m_branch_truthJet_back).c_str(), &(ED_truthJet->mjj));

      tree->SetBranchStatus  ((m_branch_truthJet_front + "m23" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "m23" + m_branch_truthJet_back).c_str(), &(ED_truthJet->m23));

      if(CS->doyStar){
      tree->SetBranchStatus  ((m_branch_truthJet_front + "yStar" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "yStar" + m_branch_truthJet_back).c_str(), &(ED_truthJet->yStar));
      }
      if(CS->doDeltaPhi){
      tree->SetBranchStatus  ((m_branch_truthJet_front + "deltaPhi" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "deltaPhi" + m_branch_truthJet_back).c_str(), &(ED_truthJet->deltaPhi));
      }
      if(CS->doPTBalance){
      tree->SetBranchStatus  ((m_branch_truthJet_front + "pTBalance" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "pTBalance" + m_branch_truthJet_back).c_str(), &(ED_truthJet->pTBalance));
      }
      if(CS->doMHT){
      tree->SetBranchStatus  ((m_branch_truthJet_front + "MHT" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "MHT" + m_branch_truthJet_back).c_str(), &(ED_truthJet->MHT));
      }
      if(CS->doMHTPhi){
      tree->SetBranchStatus  ((m_branch_truthJet_front + "MHTPhi" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "MHTPhi" + m_branch_truthJet_back).c_str(), &(ED_truthJet->MHTPhi));
      }
      if(CS->doEMFrac){
      tree->SetBranchStatus  ((m_branch_truthJet_front + "EMFrac" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "EMFrac" + m_branch_truthJet_back).c_str(), &(ED_truthJet->EMFrac));
      }
      if(CS->doHECFrac){
      tree->SetBranchStatus  ((m_branch_truthJet_front + "HECFrac" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "HECFrac" + m_branch_truthJet_back).c_str(), &(ED_truthJet->HECFrac));
      }
      if(CS->doFracSamplingMax){
      tree->SetBranchStatus  ((m_branch_truthJet_front + "FracSamplingMax" + m_branch_truthJet_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet_front + "FracSamplingMax" + m_branch_truthJet_back).c_str(), &(ED_truthJet->FracSamplingMax));
      }

  }

  // L1 jets
  if ((CS->doTurnOns)&&(CS->useEmulation)){

      if(m_debug) Info("InputHandler()", "Setting L1 branches");

      tree->SetBranchStatus  ((m_branch_LVL1JetROI_front + "eta" + m_branch_LVL1JetROI_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_LVL1JetROI_front + "eta" + m_branch_LVL1JetROI_back).c_str(), &(L1D->eta));

      tree->SetBranchStatus  ((m_branch_LVL1JetROI_front + "phi" + m_branch_LVL1JetROI_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_LVL1JetROI_front + "phi" + m_branch_LVL1JetROI_back).c_str(), &(L1D->phi));

      tree->SetBranchStatus  ((m_branch_LVL1JetROI_front + "et8x8" + m_branch_LVL1JetROI_back).c_str(), 1);
      tree->SetBranchAddress ((m_branch_LVL1JetROI_front + "et8x8" + m_branch_LVL1JetROI_back).c_str(), &(L1D->Et));

  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode InputHandler :: initialize ()
{
  if(m_debug) std::cout << "Calling InputHandler :: initialize" << std::endl;
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.


  // Initialise GRL
  if ((CS->isData) && (m_applyGRL)){
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    std::vector<std::string> vecStringGRL;
    m_GRLxml = gSystem->ExpandPathName( m_GRLxml.c_str() );
    vecStringGRL.push_back(m_GRLxml);
    RETURN_CHECK("BasicEventSelection::initialize()", m_grl->setProperty( "GoodRunsListVec", vecStringGRL), "");
    RETURN_CHECK("BasicEventSelection::initialize()", m_grl->setProperty("PassThrough", false), "");
    RETURN_CHECK("BasicEventSelection::initialize()", m_grl->initialize(), "");
  }

  Info("initialize()", "Succesfully initialized! \n");
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode InputHandler :: execute ()
{
  if(m_debugInExecute) Info("InputHandler()", "Calling InputHandler :: execute");

  // Here you do everything that needs to be done on every single
  // event, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // EVENTCOUNTER
  if (CS->doOnlyThisNumberOfEvents > -1){
    if(eventCounter >= CS->doOnlyThisNumberOfEvents){
        eventCounter++;
        return EL::StatusCode::SUCCESS;
    }
  }

  if (m_coutPassedTriggers){
      cout << "=== Passed Triggers of event: " << eventCounter << endl;
      for (int n=0; n < TD->event_passedTriggers->size(); n++){
          cout << TD->event_passedTriggers->at(n) << endl;
      }
  }
  /*  // Crosscheck Timing!!!! TODO
  cout << "Timing of event: " << eventCounter << endl;
  for (unsigned int n=0; n < ED_jet->Timing->size(); n++){
      cout << "ED_jet->Timing->at(" << n << "): " << ED_jet->Timing->at(n) << endl;
  }
  for (unsigned int n=0; n < ED_trigJet->Timing->size(); n++){
      cout << "ED_trigJet->Timing->at(" << n << "): " << ED_trigJet->Timing->at(n) << endl;
      }*/


  //Get entries out of the tree
  wk()->tree()->GetEntry (wk()->treeEntry());

  //Set E and pt vector for L1Data
  if ((CS->doTurnOns)&&(CS->useEmulation)) L1D->FillEAndPTVector();

  //For MC only: Normalise weight by dividing through numberOfEntries
  if (!CS->isData){
      m_weight = m_weight/ (double) numberOfEntries;
  }


  // --- check if jet vector has more than 1 jets otherwise skip this event
  // check offline jets
  if ((CS->doOfflineJetKinematics)||(CS->doTriggerTruthResponse)||(CS->doTriggerOfflineResponse)|| (CS->doMjjResponseOffVsTruth) || (CS->doMjjResponseTrigVsTruth) || (CS->doTurnOns)){
      if(ED_jet->eta->size() < 2) return EL::StatusCode::SUCCESS;
  }
  // check trigger jets
  if ((CS->doTriggerJetKinematics)||(CS->doTriggerTruthResponse)||(CS->doTriggerOfflineResponse)|| (CS->doMjjResponseTrigVsOff) || (CS->doMjjResponseTrigVsTruth) || ((CS->doTurnOns)&&(CS->useEmulation))){
      if(ED_trigJet->eta->size() < 2) return EL::StatusCode::SUCCESS;
  }
  // check truth jets
  if ((!CS->isData) &&
      ((CS->doTruthJetKinematics)||(CS->doOfflineTruthResponse)|| (CS->doMjjResponseOffVsTruth) || (CS->doMjjResponseTrigVsTruth) || (CS->doTriggerTruthResponse))){
      if(ED_truthJet->eta->size() < 2) return EL::StatusCode::SUCCESS;
  }

  // GRL
  if ((CS->isData) && (m_applyGRL)){

    if (!m_grl->passRunLB( m_runNumber, m_lumiBlock ) ) {
        if(m_debugInExecute) cout << "GRL:: Fail Event " << endl;
        return EL::StatusCode::SUCCESS; // go to next event
    }else{
        if(m_debugInExecute) cout << "GRL:: Pass Event " << endl;
    }

  }

  if (m_debugInExecute){
      cout << "eventCounter: " << eventCounter << endl;
  }


  /* ======= CUTS ======= */

  // Reset Cut Class
  cutH->Reset(TD->config_passedTriggers);
  cutH->Reset(TD->probe_passedTrigger);
  cutH->Reset(TD->ref_passedTrigger);

  // CHECK LArERROR
  if ((CS->isData) && (CS->checkLArError)) cutH->AddCut(!m_LArError); // if argument is false, event does NOT pass

  if (CS->doCuts){
      // Add simple cuts on observables
      //possible options:
      //"max": cutValue is maximal value
      //"min": cutValue is minimal value
      //"absmax": cutValue is maximal value of |observable|
      //"absmin": cutValue is minimal value of |observable|
      // cuts in GeV!

      // leading pt cut
      //cutH->AddCut(ED_trigJet->pt->at(0), 30.0, "min");

      // subleading pt cut
      //cutH->AddCut(ED_trigJet->pt->at(1), 85.0, "min");

      // eta cut on leading and subleading jet
      //cutH->AddCut(ED_jet->eta->at(0), 2.8, "absmax"); //OFFLINE JETS
      //cutH->AddCut(ED_trigJet->eta->at(1), 2.8, "absmax");

      // yStar cut
      //cutH->AddCut(ED_trigJet->yStar, 0.6, "absmax");

      // mjj cut
      //cutH->AddCut(ED_trigJet->mjj, 500.0, "min");
      //cutH->AddCut(ED_trigJet->mjj, 900.0, "max");

      // cutString method
      cutH->UseCutStringMethod(ED_jet, ED_trigJet, ED_truthJet);

  }

  // Discard if event did not pass
  if (!cutH->passedCuts()) {return EL::StatusCode::SUCCESS;}

  // Investigate which triggers this event passed and right the result out in a bool vector
  // TD->config_triggerName : string vector of wanted triggers
  // TD->event_passedTriggers : string vector of triggers that this event passed
  // TD->config_passedTriggers: bool that gives the status of wanted triggers for this event

  cutH->SearchTrigger(TD->config_triggerName, TD->event_passedTriggers, TD->config_passedTriggers);

  // Fill Turnon relied booleans
  if (CS->doTurnOns){

      // probe
      cutH->SearchTrigger(TD->probe_triggerName, TD->event_passedTriggers, TD->probe_passedTrigger);

      // ref
      cutH->SearchTrigger(TD->ref_triggerName, TD->event_passedTriggers, TD->ref_passedTrigger);

      if (CS->useTriggerBeforePraescale){

          //probe
          cutH->SearchPassBits(TD->probe_triggerName,TD->probe_passBits, TD->event_isPassBitsNames, TD->event_isPassBits);

          //ref
          cutH->SearchPassBits(TD->ref_triggerName,TD->ref_passBits, TD->event_isPassBitsNames, TD->event_isPassBits);
      }

  }

  //Apply trigger specific cuts
  /*  if (CS->doCuts){
      cutH->ApplyTriggerSpecificCut("HLT_j360", TD->config_triggerName, TD->config_passedTriggers, ED_trigJet->pt->at(0), 410, "min");
      cutH->ApplyTriggerSpecificCut("HLT_j260", TD->config_triggerName, TD->config_passedTriggers, ED_trigJet->pt->at(0), 280, "min");
      cutH->ApplyTriggerSpecificCut("HLT_j60", TD->config_triggerName, TD->config_passedTriggers, ED_trigJet->pt->at(0), 85, "min");
      cutH->ApplyTriggerSpecificCut("HLT_j85", TD->config_triggerName, TD->config_passedTriggers, ED_trigJet->pt->at(0), 110, "min");
      } */

  if (m_debugInExecute){
      cout << "Event passed cuts!" << endl;
      cout << "The triggers that are wanted:..." << TD->config_triggerName.size() << endl;
      for (int n=0; n < TD->config_triggerName.size(); n++){
        cout << TD->config_triggerName.at(n) << " passes?: " << TD->config_passedTriggers.at(n) << endl;
      }
  }


  /* ======= Matching ======= */

  // 1.1. Matching of trigger and offline jets

  // Initialise matchingIndex
  std::vector<int> matchingIndexList_TriggvsOffl;
  std::vector<int> matchingIndexList_OfflvsTruth;
  std::vector<int> matchingIndexList_TriggvsTruth;

  std::vector<double> DeltaRList_TriggvsOffl;
  std::vector<double> DeltaRList_OfflvsTruth;
  std::vector<double> DeltaRList_TriggvsTruth;

  // Just do matching if it is really needed, i.e. if we want to see Response plots
  if (CS->doTriggerOfflineResponse || CS->doOfflineTruthResponse || CS->doTriggerTruthResponse || CS->doMjjResponseOffVsTruth || CS->doMjjResponseTrigVsOff || CS->doMjjResponseTrigVsTruth){

    if (CS->doMatching){
        for (unsigned int i = 0; i < ED_trigJet->pt->size(); i++){
            matchingIndexList_TriggvsOffl.push_back(-1);
            matchingIndexList_TriggvsTruth.push_back(-1);

            DeltaRList_TriggvsOffl.push_back(-1);
            DeltaRList_TriggvsTruth.push_back(-1);
        }

        for (unsigned int i = 0; i < ED_jet->pt->size(); i++){
            matchingIndexList_OfflvsTruth.push_back(-1);
            DeltaRList_OfflvsTruth.push_back(-1);
        }

        // Matching of offline and trigger jets
        if ((CS->doTriggerOfflineResponse) || (CS->doMjjResponseTrigVsOff)) {
            anaHandler->findBestMatching(ED_trigJet, ED_jet, matchingIndexList_TriggvsOffl, DeltaRList_TriggvsOffl, m_DeltaRMax);
        }

        // Matching of offline and truth jets
        if ((CS->doOfflineTruthResponse) || (CS->doMjjResponseOffVsTruth)) {
            anaHandler->findBestMatching(ED_jet, ED_truthJet, matchingIndexList_OfflvsTruth, DeltaRList_OfflvsTruth, m_DeltaRMax);
        }

        // Matching of trigger and truht jets
        if ((CS->doTriggerTruthResponse) || (CS->doMjjResponseTrigVsTruth)){
            anaHandler->findBestMatching(ED_trigJet, ED_truthJet, matchingIndexList_TriggvsTruth, DeltaRList_TriggvsTruth, m_DeltaRMax);
        }
    }
    else{
        // if matchingIndex is still used, give it reasonable entries
        int counter = 0;
        for (unsigned int i = 0; i < ED_trigJet->pt->size(); i++){
          matchingIndexList_TriggvsOffl.push_back(counter);
          matchingIndexList_TriggvsTruth.push_back(counter);
          counter ++;
          DeltaRList_TriggvsOffl.push_back(-1);
          DeltaRList_TriggvsTruth.push_back(-1);
        }


        for (unsigned int i = 0; i < ED_jet->pt->size(); i++){
            matchingIndexList_OfflvsTruth.push_back(counter);
            counter++;
            DeltaRList_OfflvsTruth.push_back(-1);
        }
    }
  }

  /* ======= CHECK JET BY JET => isGood ======= */

  if (m_debugInExecute) cout << "before isGood.." << endl;

  // OFFLINE JETS
  if ((CS->doOfflineJetKinematics)||(CS->doOfflineTruthResponse)||(CS->doTriggerOfflineResponse)){
      for (unsigned int i = 0; i < ED_jet->pt->size(); i++) {
        // first set true
        ED_jet->isGood[i] = true;
        if ((CS->isData) && (CS->doCleaning)) {
	      if (ED_jet->passedCleaning->at(i) == 0) ED_jet->isGood[i] = false;
        }
      }
  }

  // TRIGGER JETS
  if ((CS->doTriggerJetKinematics)||(CS->doTriggerTruthResponse)||(CS->doTriggerOfflineResponse)){
      for (unsigned int i = 0; i < ED_trigJet->pt->size(); i++) {
        // first set true
        ED_trigJet->isGood[i] = true;
        if ((CS->isData) && (CS->doCleaning)) {
	      if (ED_trigJet->passedCleaning->at(i) == 0) ED_trigJet->isGood[i] = false;
        }
      }
  }

  // TRUTH JETS (No cleaning performed here, but isGood vector is needed for filling)
  if ((CS->doTruthJetKinematics)||(CS->doOfflineTruthResponse)||(CS->doTriggerTruthResponse)){
      for (unsigned int i = 0; i < ED_truthJet->pt->size(); i++) {
        // first set true
        ED_truthJet->isGood[i] = true;
        }
  }

  if (m_debugInExecute) cout << "after isGood.." << endl;

  /* ======= Fill Kinematics ======= */

  for (unsigned int n=0; n < TD->config_passedTriggers.size(); n++){
      if (m_debugInExecute) cout << "trigger: " << TD->config_triggerName.at(n) << " ::: bool value (passed this event?): "<< TD->config_passedTriggers.at(n) << endl;
      if (!TD->config_passedTriggers.at(n)) continue;

      // Fill triggerHistoPacks
      m_map_triggerhistos[TD->config_triggerName.at(n)]->FillAll(ED_jet, ED_trigJet, ED_truthJet, matchingIndexList_TriggvsOffl, matchingIndexList_OfflvsTruth, matchingIndexList_TriggvsTruth, m_weight);
  }

  /* ======= Fill Turnons ======= */
  if (CS->doTurnOns){
      trigEffic->FillAll(TD, ED_jet, ED_trigJet, ED_truthJet, L1D, m_weight, CS);
  }


  eventCounter++;
  TD->UpdateCounting();
  if (m_debug) cout << "=== End of Execute ===" << endl;
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode InputHandler :: postExecute ()
{
  //if(m_debug)  std::cout << "EDDI: Calling InputHandler :: postExecute" << std::endl;
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode InputHandler :: finalize ()
{
  if(m_debug) Info("InputHandler()", "Calling InputHandler :: finalize");
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  // Creating Turn-On histograms
  if (CS->doTurnOns){
    trigEffic->DivideEfficiencyPlots(TD, CS, wk());
  }

  cout << "================== Results ====================" << endl;

  std::cout << "Number of events in dataset: " << eventCounter << std::endl;

  TD->CoutCounting();

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode InputHandler :: histFinalize ()
{
  if(m_debug) Info("InputHandler()", "Calling InputHandler :: histFinalize");
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.

  out_trigger->close();
  out_logfile->close();
  logwr->CloseLogfile();

  return EL::StatusCode::SUCCESS;
}

