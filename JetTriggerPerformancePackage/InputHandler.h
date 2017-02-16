/*--------------------------------------------------------------------
InputHandler.h

The InputHandler is JTPP's main function where all wanted histogram classes
are declared and their methods are used.

created by Edgar Kellermann (edgar.kellermann@cern.ch)
--------------------------------------------------------------------*/

#ifndef JetTriggerPerformancePackage_InputHandler_H
#define JetTriggerPerformancePackage_InputHandler_H

#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>
#include <EventLoop/Worker.h>

// rootcore includes
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

//algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// ROOT include(s):
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TLorentzVector.h"

#include <sstream>
#include <vector>

// own classes
#include <JetTriggerPerformancePackage/TriggerHistoPack.h>
#include <JetTriggerPerformancePackage/ResponseMatrix.h>
#include <JetTriggerPerformancePackage/KinematicMatrix.h>
#include <JetTriggerPerformancePackage/TriggerEfficiencyMatrix.h>
#include <JetTriggerPerformancePackage/ToolsJTPP.h>
#include <JetTriggerPerformancePackage/AnalysisHandler.h>
#include <JetTriggerPerformancePackage/LogWriter.h>
#include <JetTriggerPerformancePackage/EventData.h>
#include <JetTriggerPerformancePackage/L1Data.h>
#include <JetTriggerPerformancePackage/TriggerData.h>
#include <JetTriggerPerformancePackage/ConfigStatus.h>
#include <JetTriggerPerformancePackage/CutHandler.h>


using namespace std;



class InputHandler : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  public:

    // GRL
    bool m_applyGRL;
    string m_GRLxml;

    // Debug
    bool m_debug;
    bool m_debugInExecute;
    bool m_coutPassedTriggers;

    // branch strings
    string m_branch_offlineJet;
    string m_branch_triggerJet;
    string m_branch_truthJet;
    string m_branch_jet_truth;
    string m_branch_LVL1JetROI;

    // specify triggers of interest
    string m_TriggerName;
    vector<string> m_confTriggers;

    bool m_isData;
    bool m_checkLArError;

    bool m_doOfflineJetKinematics;
    bool m_doTriggerJetKinematics;
    bool m_doTruthJetKinematics;

    bool m_doLeadSubleadThirdJetKinematics;
    bool m_doNthJetKinematics;
    bool m_doKinematicsInBinsOfEta;
    bool m_doEtainBinsOfPt;

    bool m_doMjj;
    bool m_doM23;
    bool m_doyStar;
    bool m_doDeltaPhi;
    bool m_doPTBalance;
    bool m_doHT;
    bool m_doMHT;
    bool m_doMHTPhi;
    bool m_doEMFrac;
    bool m_doHECFrac;
    bool m_doFracSamplingMax;

    bool m_calculateMjj;
    bool m_calculateM23;
    bool m_calculateYStar;

    string m_Kinematic_E_Binning;
    string m_Kinematic_pt_Binning;
    string m_Kinematic_phi_Binning;
    string m_Kinematic_eta_Binning;

    string m_Kinematic_mjj_Binning;
    string m_Kinematic_m23_Binning;
    string m_Kinematic_yStar_Binning;
    string m_Kinematic_deltaPhi_Binning;
    string m_Kinematic_pTBalance_Binning;
    string m_Kinematic_HT_Binning;
    string m_Kinematic_MHT_Binning;
    string m_Kinematic_MHTPhi_Binning;
    string m_Kinematic_EMFrac_Binning;
    string m_Kinematic_HECFrac_Binning;
    string m_Kinematic_FracSamplingMax_Binning;

    bool m_doMatching;
    float m_DeltaRMax;

    bool m_doOfflineTruthResponse;
    bool m_doTriggerTruthResponse;
    bool m_doTriggerOfflineResponse;

    bool m_doMjjResponseTrigVsOff;
    bool m_doMjjResponseOffVsTruth;
    bool m_doMjjResponseTrigVsTruth;

    int m_Response_BinNumbers;
    float m_Response_MinBin;
    float m_Response_MaxBin;

    string m_PtResponse_ptBinning;
    string m_PtResponse_etaBinning;
    string m_MjjResponse_mjjBinning;
    string m_MjjResponse_etaBinning;

    bool m_doTurnOns;
    string m_TurnOnName;
    bool m_useTriggerDecisionTool;
    bool m_useEmulation;
    bool m_useTriggerBeforePraescale;

    int   m_PtTurnon_BinNumbers;
    float m_PtTurnon_MinBin;
    float m_PtTurnon_MaxBin;

    int   m_HtTurnon_BinNumbers;
    float m_HtTurnon_MinBin;
    float m_HtTurnon_MaxBin;

    int   m_MjjTurnon_BinNumbers;
    float m_MjjTurnon_MinBin;
    float m_MjjTurnon_MaxBin;

    // Turnon Event Selection
    float m_HLT_cen_etaMin;
    float m_HLT_cen_etaMax;
    float m_HLT_fwd_etaMin;
    float m_HLT_fwd_etaMax;
    float m_L1_cen_etaMin;
    float m_L1_cen_etaMax;
    float m_L1_fwd_etaMin;
    float m_L1_fwd_etaMax;
    float m_TurnonCut_Timing;

    std::string m_cutStringTurnons;

    // Kinematic Cuts
    bool m_doCuts;
    std::string m_cutStringKinematics;
    bool m_doCleaning;

    int m_doOnlyThisNumberOfEvents;

    std::string m_name;
    float m_mcEventWeight;  //!

    int* m_testArray;


  private:

    // Branch Names
    string m_branch_offlineJet_front;
    string m_branch_triggerJet_front;
    string m_branch_truthJet_front;
    string m_branch_jet_truth_front;
    string m_branch_LVL1JetROI_front;

    string m_branch_offlineJet_back;
    string m_branch_triggerJet_back;
    string m_branch_truthJet_back;
    string m_branch_jet_truth_back;
    string m_branch_LVL1JetROI_back;

    int numberOfEntries;

    GoodRunsListSelectionTool*   m_grl;       //!

    int eventCounter; // 'It's the final countdown!'

    int m_runNumber; //!
    Long64_t m_eventNumber; //!
    int m_lumiBlock; //!
    int m_mcEventNumber; //!
    int m_ChannelNumber; //!

    // Data specific variables
    bool m_LArError; //!

    // variables for event observables
    float m_weight; //!

    // ConfigStatus
    ConfigStatus* CS; //!

    //TriggerHistoPack
    TriggerHistoPack* triggerhistos;
    std::map<std::string, TriggerHistoPack*> m_map_triggerhistos;

    //TriggerEfficiencyMatrix
    TriggerEfficiencyMatrix* trigEffic;

    // EventData
    EventData* ED_jet; //!
    EventData* ED_trigJet; //!
    EventData* ED_truthJet; //!
    EventData* ED_jet_truth; //!
    L1Data* L1D; //!

    // TriggerData
    TriggerData* TD; //!

    // CutHandler
    CutHandler* cutH;

    // AnalysisHandler
    AnalysisHandler* anaHandler; //!

    // LogWriter
    LogWriter* logwr;

    // output of logfiles for for JTPPPlotter.py
    ofstream* out_trigger; //!
    ofstream* out_turnOns; //!
    ofstream* out_logfile; //!
    //    ofstream out_leadSubleadThird;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  // this is a standard constructor
  InputHandler();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // these are the functions not inherited from Algorithm
  virtual EL::StatusCode configure ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(InputHandler, 1);
};

#endif
