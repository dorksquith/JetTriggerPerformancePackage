/*--------------------------------------------------------------------
ConfigStatus.h

ConfigStatus is meant to be a container for all settings  that were made in the config file.

created by Edgar Kellermann (edgar.kellermann@cern.ch)
--------------------------------------------------------------------*/

#ifndef JetTriggerPerformancePackage_CONFIGSTATUS_H
#define JetTriggerPerformancePackage_CONFIGSTATUS_H

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "TLorentzVector.h"


//Own classes
#include <JetTriggerPerformancePackage/BaseJTPP.h>

class ConfigStatus: public BaseJTPP
{
 public:
    ConfigStatus(bool a_debug,
		 bool a_debugInExecute,
		 std::string a_TriggerName,
		 bool a_isData, bool a_checkLArError,
		 bool a_doOfflineJetKinematics,
		 bool a_doTriggerJetKinematics,
		 bool a_doTruthJetKinematics,
		 bool a_doLeadSubleadThirdJetKinematics,
		 bool a_doNthJetKinematics,
		 bool a_doKinematicsInBinsOfEta,
		 bool a_doEtainBinsOfPt,
         bool a_doMjj,
         bool a_doM23,
		 bool a_doyStar,
		 bool a_doDeltaPhi,
		 bool a_doPTBalance,
         bool a_doHT,
		 bool a_doMHT,
		 bool a_doMHTPhi,
		 bool a_doEMFrac,
		 bool a_doHECFrac,
		 bool a_doFracSamplingMax,
         bool a_calculateMjj,
         bool a_calculateM23,
         bool a_calculateYStar,
         std::string a_Kinematic_E_Binning,
         std::string a_Kinematic_pt_Binning,
         std::string a_Kinematic_phi_Binning,
         std::string a_Kinematic_eta_Binning,
         std::string a_Kinematic_mjj_Binning,
         std::string a_Kinematic_m23_Binning,
         std::string a_Kinematic_yStar_Binning,
         std::string a_Kinematic_deltaPhi_Binning,
         std::string a_Kinematic_pTBalance_Binning,
         std::string a_Kinematic_HT_Binning,
         std::string a_Kinematic_MHT_Binning,
         std::string a_Kinematic_MHTPhi_Binning,
         std::string a_Kinematic_EMFrac_Binning,
         std::string a_Kinematic_HECFrac_Binning,
         std::string a_Kinematic_FracSamplingMax_Binning,
         bool a_doMatching,
         bool a_doOfflineTruthResponse,
		 bool a_doTriggerTruthResponse,
		 bool a_doTriggerOfflineResponse,
		 bool a_doMjjResponseTrigVsOff,
		 bool a_doMjjResponseOffVsTruth,
		 bool a_doMjjResponseTrigVsTruth,
         int a_Response_BinNumbers,
         float a_Response_MinBin,
         float a_Response_MaxBin,
         std::string a_PtResponse_ptBinning,
         std::string a_PtResponse_etaBinning,
         std::string a_MjjResponse_mjjBinning,
         std::string a_MjjResponse_etaBinning,
		 bool a_doTurnOns,
         std::string a_TurnOnName,
		 bool a_useTriggerDecisionTool,
		 bool a_useEmulation,
		 bool a_useTriggerBeforePraescale,
         int   a_PtTurnon_BinNumbers,
         float a_PtTurnon_MinBin,
         float a_PtTurnon_MaxBin,
         int   a_HtTurnon_BinNumbers,
         float a_HtTurnon_MinBin,
         float a_HtTurnon_MaxBin,
         int   a_MjjTurnon_BinNumbers,
         float a_MjjTurnon_MinBin,
         float a_MjjTurnon_MaxBin,
         const float a_HLT_cen_etaMin,
         const float a_HLT_cen_etaMax,
         const float a_HLT_fwd_etaMin,
         const float a_HLT_fwd_etaMax,
         const float a_L1_cen_etaMin,
         const float a_L1_cen_etaMax,
         const float a_L1_fwd_etaMin,
         const float a_L1_fwd_etaMax,
         const float a_TurnonCut_Timing,
		 bool a_doCuts,
         std::string a_cutStringKinematics,
         std::string a_cutStringTurnons,
		 bool a_doCleaning,
		 int a_doOnlyThisNumberOfEvents); // "a_" for argument
    virtual ~ConfigStatus();

    // Debug
    const bool debug;
    const bool debugInExecute;

    // Trigger Selection
    const std::string TriggerName;

    // Data or MC
    const bool isData;

    //For Data Samples
    const bool checkLArError;

    //Kinematics
    const bool doOfflineJetKinematics;
    const bool doTriggerJetKinematics;
    const bool doTruthJetKinematics;

    //More setups for kinematics
    const bool doLeadSubleadThirdJetKinematics;
    const bool doNthJetKinematics;
    const bool doKinematicsInBinsOfEta;
    const bool doEtainBinsOfPt;

    //More kinematics
    const bool doMjj;
    const bool doM23;
    const bool doyStar;
    const bool doDeltaPhi;
    const bool doPTBalance;
    const bool doHT;
    const bool doMHT;
    const bool doMHTPhi;
    const bool doEMFrac;
    const bool doHECFrac;
    const bool doFracSamplingMax;

    //Calculate observables
    const bool calculateMjj;
    const bool calculateM23;
    const bool calculateYStar;

    // Kinematic Binning
    const std::string Kinematic_E_Binning;
    const std::string Kinematic_pt_Binning;
    const std::string Kinematic_phi_Binning;
    const std::string Kinematic_eta_Binning;
    const std::string Kinematic_mjj_Binning;
    const std::string Kinematic_m23_Binning;
    const std::string Kinematic_yStar_Binning;
    const std::string Kinematic_deltaPhi_Binning;
    const std::string Kinematic_pTBalance_Binning;
    const std::string Kinematic_HT_Binning;
    const std::string Kinematic_MHT_Binning;
    const std::string Kinematic_MHTPhi_Binning;
    const std::string Kinematic_EMFrac_Binning;
    const std::string Kinematic_HECFrac_Binning;
    const std::string Kinematic_FracSamplingMax_Binning;

    // Matching
    const bool doMatching;

    //ptResponse
    const bool doOfflineTruthResponse;
    const bool doTriggerTruthResponse;
    const bool doTriggerOfflineResponse;

    //mjjResponse
    const bool doMjjResponseTrigVsOff;
    const bool doMjjResponseOffVsTruth;
    const bool doMjjResponseTrigVsTruth;

    //Response Binning
    const int   Response_BinNumbers;
    const float Response_MinBin;
    const float Response_MaxBin;

    //Pt and Eta Response Binnings
    const std::string PtResponse_ptBinning;
    const std::string PtResponse_etaBinning;
    const std::string MjjResponse_mjjBinning;
    const std::string MjjResponse_etaBinning;

    //turn-ons
    const bool doTurnOns;
    const std::string TurnOnName;
    const bool useTriggerDecisionTool;
    const bool useEmulation;
    const bool useTriggerBeforePraescale;
    const int   PtTurnon_BinNumbers;
    const float PtTurnon_MinBin;
    const float PtTurnon_MaxBin;
    const int   HtTurnon_BinNumbers;
    const float HtTurnon_MinBin;
    const float HtTurnon_MaxBin;
    const int   MjjTurnon_BinNumbers;
    const float MjjTurnon_MinBin;
    const float MjjTurnon_MaxBin;
    const float HLT_cen_etaMin;
    const float HLT_cen_etaMax;
    const float HLT_fwd_etaMin;
    const float HLT_fwd_etaMax;
    const float L1_cen_etaMin;
    const float L1_cen_etaMax;
    const float L1_fwd_etaMin;
    const float L1_fwd_etaMax;
    const float TurnonCut_Timing;

    //cuts and cleaning
    const bool doCuts;
    const std::string cutStringKinematics;
    const std::string cutStringTurnons;
    const bool doCleaning;
    const int doOnlyThisNumberOfEvents;

};

#endif // JetTriggerPerformancePackage_CONFIGSTATUS_H
