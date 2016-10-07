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
    m_branch_offlineJet("jet"),
    m_branch_triggerJet("trigJet"),
    m_branch_truthJet("truthJet"),
    m_branch_jet_truth("jet_truth"),
    m_branch_LVL1JetROI("LVL1JetROI"),
    m_TriggerName("STUDYALL"),
    m_debug(false),
    m_debugInExecute(false),
    m_isData(false),
    m_checkLArError(false),
    m_doOfflineJetKinematics(false),
    m_doTriggerJetKinematics(false),
    m_doTruthJetKinematics(false),
    m_doLeadSubleadThirdJetKinematics(false),
    m_doNthJetKinematics(false),
    m_doKinematicsInBinsOfEta(false),
    m_doEtainBinsOfPt(false),
    m_doyStar(false),
    m_doDeltaPhi(false),
    m_doPTBalance(false),
    m_doMHT(false),
    m_doMHTPhi(false),
    m_doEMFrac(false),
    m_doHECFrac(false),
    m_doFracSamplingMax(false),
    m_doOfflineTruthResponse(false),
    m_doTriggerTruthResponse(false),
    m_doTriggerOfflineResponse(false),
    m_doMjjResponseTrigVsOff(false),
    m_doMjjResponseOffVsTruth(false),
    m_doMjjResponseTrigVsTruth(false),
    m_doTurnOns(false),
    m_TurnOnName("STUDYALL/STUDYALL"),
    m_useTriggerDecisionTool(false),
    m_useEmulation(false),
    m_useTriggerBeforePraescale(false),
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
    ED_jet_truth(nullptr),
    L1D(nullptr),
    anaHandler(nullptr),
    logwr(nullptr),
    DeltaRMax(0.4),
    out_trigger(nullptr),
    out_turnOns(nullptr),
    out_logfile(nullptr)
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
                        m_doyStar,
                        m_doDeltaPhi,
                        m_doPTBalance,
                        m_doMHT,
                        m_doMHTPhi,
                        m_doEMFrac,
                        m_doHECFrac,
                        m_doFracSamplingMax,
                        m_doOfflineTruthResponse,
                        m_doTriggerTruthResponse,
                        m_doTriggerOfflineResponse,
                        m_doMjjResponseTrigVsOff,
                        m_doMjjResponseOffVsTruth,
                        m_doMjjResponseTrigVsTruth,
                        m_doTurnOns,
                        m_TurnOnName,
                        m_useTriggerDecisionTool,
                        m_useEmulation,
                        m_useTriggerBeforePraescale,
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
  std::cout << "\n==== Kinmatics Cut Selection ====" << std::endl;
  cutH = new CutHandler(CS->cutStringKinematics);

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
  ED_jet = new EventData("jet");
  ED_trigJet = new EventData("trigJet");
  ED_truthJet = new EventData("truthJet");
  ED_jet_truth = new EventData("jet_truth");
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

  tree->SetBranchStatus  ("runNumber",    1);
  tree->SetBranchAddress ("runNumber",    &m_runNumber);

  tree->SetBranchStatus  ("eventNumber",    1);
  tree->SetBranchAddress ("eventNumber",    &m_eventNumber);

  tree->SetBranchStatus  ("weight",    1);
  tree->SetBranchAddress ("weight",    &m_weight);


  // Data specific observables
  if (CS->checkLArError){
      tree->SetBranchStatus  ("LArError", 1);
      tree->SetBranchAddress ("LArError", &m_LArError);
  }

  // Offline jets
  if ((CS->doOfflineJetKinematics)||(CS->doOfflineTruthResponse)||(CS->doTriggerOfflineResponse)||(CS->doTurnOns)){

      if(m_debug) Info("InputHandler()", "Setting offline branches");

      tree->SetBranchStatus  ((m_branch_offlineJet + "_pt").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_pt").c_str(), &(ED_jet->pt));

      tree->SetBranchStatus  ((m_branch_offlineJet + "_eta").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_eta").c_str(), &(ED_jet->eta));

      tree->SetBranchStatus  ((m_branch_offlineJet + "_phi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_phi").c_str(), &(ED_jet->phi));

      tree->SetBranchStatus  ((m_branch_offlineJet + "_E").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_E").c_str(), &(ED_jet->E));

      tree->SetBranchStatus  ((m_branch_offlineJet + "_mjj").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_mjj").c_str(), &(ED_jet->mjj));

      tree->SetBranchStatus  ((m_branch_offlineJet + "_m23").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_m23").c_str(), &(ED_jet->m23));

      if (CS->isData){
      tree->SetBranchStatus  ((m_branch_offlineJet + "_clean_passLooseBad").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_clean_passLooseBad").c_str(), &(ED_jet->passedCleaning));
      }

      if(CS->doyStar){
      tree->SetBranchStatus  ((m_branch_offlineJet + "_yStar").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_yStar").c_str(), &(ED_jet->yStar));
      }
      if(CS->doDeltaPhi){
      tree->SetBranchStatus  ((m_branch_offlineJet + "_deltaPhi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_deltaPhi").c_str(), &(ED_jet->deltaPhi));
      }
      if(CS->doPTBalance){
      tree->SetBranchStatus  ((m_branch_offlineJet + "_pTBalance").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_pTBalance").c_str(), &(ED_jet->pTBalance));
      }
      if(CS->doMHT){
      tree->SetBranchStatus  ((m_branch_offlineJet + "_MHT").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_MHT").c_str(), &(ED_jet->MHT));
      }
      if(CS->doMHTPhi){
      tree->SetBranchStatus  ((m_branch_offlineJet + "_MHTPhi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_MHTPhi").c_str(), &(ED_jet->MHTPhi));
      }
      if(CS->doEMFrac){
      tree->SetBranchStatus  ((m_branch_offlineJet + "_EMFrac").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_EMFrac").c_str(), &(ED_jet->EMFrac));
      }
      if(CS->doHECFrac){
      tree->SetBranchStatus  ((m_branch_offlineJet + "_HECFrac").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_HECFrac").c_str(), &(ED_jet->HECFrac));
      }
      if(CS->doFracSamplingMax){
      tree->SetBranchStatus  ((m_branch_offlineJet + "_FracSamplingMax").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_FracSamplingMax").c_str(), &(ED_jet->FracSamplingMax));
      }
      if (CS->doTurnOns){
      tree->SetBranchStatus  ((m_branch_offlineJet + "_Timing").c_str(), 1);
      tree->SetBranchAddress ((m_branch_offlineJet + "_Timing").c_str(), &(ED_jet->Timing));
      }
  }

  // Trigg jets
  if ((CS->doTriggerJetKinematics)||(CS->doTriggerTruthResponse)||(CS->doTriggerOfflineResponse)||((CS->doTurnOns)&&(CS->useEmulation))){

      if(m_debug) Info("InputHandler()", "Setting trigger branches");

      tree->SetBranchStatus  ((m_branch_triggerJet + "_pt").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_pt").c_str(), &(ED_trigJet->pt));

      tree->SetBranchStatus  ((m_branch_triggerJet + "_eta").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_eta").c_str(), &(ED_trigJet->eta));

      tree->SetBranchStatus  ((m_branch_triggerJet + "_phi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_phi").c_str(), &(ED_trigJet->phi));

      tree->SetBranchStatus  ((m_branch_triggerJet + "_E").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_E").c_str(), &(ED_trigJet->E));

      tree->SetBranchStatus  ((m_branch_triggerJet + "_mjj").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_mjj").c_str(), &(ED_trigJet->mjj));

      tree->SetBranchStatus  ((m_branch_triggerJet + "_m23").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_m23").c_str(), &(ED_trigJet->m23));

      if (CS->doCleaning){
      tree->SetBranchStatus  ((m_branch_triggerJet + "_clean_passLooseBad").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_clean_passLooseBad").c_str(), &(ED_trigJet->passedCleaning));
      }

      if(CS->doyStar){
      tree->SetBranchStatus  ((m_branch_triggerJet + "_yStar").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_yStar").c_str(), &(ED_trigJet->yStar));
      }
      if(CS->doDeltaPhi){
      tree->SetBranchStatus  ((m_branch_triggerJet + "_deltaPhi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_deltaPhi").c_str(), &(ED_trigJet->deltaPhi));
      }
      if(CS->doPTBalance){
      tree->SetBranchStatus  ((m_branch_triggerJet + "_pTBalance").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_pTBalance").c_str(), &(ED_trigJet->pTBalance));
      }
      if(CS->doMHT){
      tree->SetBranchStatus  ((m_branch_triggerJet + "_MHT").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_MHT").c_str(), &(ED_trigJet->MHT));
      }
      if(CS->doMHTPhi){
      tree->SetBranchStatus  ((m_branch_triggerJet + "_MHTPhi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_MHTPhi").c_str(), &(ED_trigJet->MHTPhi));
      }
      if(CS->doEMFrac){
      tree->SetBranchStatus  ((m_branch_triggerJet + "_EMFrac").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_EMFrac").c_str(), &(ED_trigJet->EMFrac));
      }
      if(CS->doHECFrac){
      tree->SetBranchStatus  ((m_branch_triggerJet + "_HECFrac").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_HECFrac").c_str(), &(ED_trigJet->HECFrac));
      }
      if(CS->doFracSamplingMax){
      tree->SetBranchStatus  ((m_branch_triggerJet + "_FracSamplingMax").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_FracSamplingMax").c_str(), &(ED_trigJet->FracSamplingMax));
      }
      if (CS->doTurnOns){
      tree->SetBranchStatus  ((m_branch_triggerJet + "_Timing").c_str(), 1);
      tree->SetBranchAddress ((m_branch_triggerJet + "_Timing").c_str(), &(ED_trigJet->Timing));
      }
  }

  // Truth jets
  if ((CS->doTruthJetKinematics)||(CS->doOfflineTruthResponse)||(CS->doTriggerTruthResponse)){

      if(m_debug) Info("InputHandler()", "Setting truth branches");

      tree->SetBranchStatus  ((m_branch_truthJet + "_pt").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_pt").c_str(), &(ED_truthJet->pt));

      tree->SetBranchStatus  ((m_branch_truthJet + "_eta").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_eta").c_str(), &(ED_truthJet->eta));

      tree->SetBranchStatus  ((m_branch_truthJet + "_phi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_phi").c_str(), &(ED_truthJet->phi));

      tree->SetBranchStatus  ((m_branch_truthJet + "_E").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_E").c_str(), &(ED_truthJet->E));

      tree->SetBranchStatus  ((m_branch_truthJet + "_mjj").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_mjj").c_str(), &(ED_truthJet->mjj));

      tree->SetBranchStatus  ((m_branch_truthJet + "_m23").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_m23").c_str(), &(ED_truthJet->m23));

      if(CS->doyStar){
      tree->SetBranchStatus  ((m_branch_truthJet + "_yStar").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_yStar").c_str(), &(ED_truthJet->yStar));
      }
      if(CS->doDeltaPhi){
      tree->SetBranchStatus  ((m_branch_truthJet + "_deltaPhi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_deltaPhi").c_str(), &(ED_truthJet->deltaPhi));
      }
      if(CS->doPTBalance){
      tree->SetBranchStatus  ((m_branch_truthJet + "_pTBalance").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_pTBalance").c_str(), &(ED_truthJet->pTBalance));
      }
      if(CS->doMHT){
      tree->SetBranchStatus  ((m_branch_truthJet + "_MHT").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_MHT").c_str(), &(ED_truthJet->MHT));
      }
      if(CS->doMHTPhi){
      tree->SetBranchStatus  ((m_branch_truthJet + "_MHTPhi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_MHTPhi").c_str(), &(ED_truthJet->MHTPhi));
      }
      if(CS->doEMFrac){
      tree->SetBranchStatus  ((m_branch_truthJet + "_EMFrac").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_EMFrac").c_str(), &(ED_truthJet->EMFrac));
      }
      if(CS->doHECFrac){
      tree->SetBranchStatus  ((m_branch_truthJet + "_HECFrac").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_HECFrac").c_str(), &(ED_truthJet->HECFrac));
      }
      if(CS->doFracSamplingMax){
      tree->SetBranchStatus  ((m_branch_truthJet + "_FracSamplingMax").c_str(), 1);
      tree->SetBranchAddress ((m_branch_truthJet + "_FracSamplingMax").c_str(), &(ED_truthJet->FracSamplingMax));
      }

  }

  // Truth jets matched to offline jets
  if ((CS->doTruthJetKinematics)||(CS->doOfflineTruthResponse)||(CS->doTriggerTruthResponse)){

      if(m_debug) Info("InputHandler()", "Setting jet_truth branches");

      tree->SetBranchStatus  ((m_branch_jet_truth + "_pt").c_str(), 1);
      tree->SetBranchAddress ((m_branch_jet_truth + "_pt").c_str(), &(ED_jet_truth->pt));

      tree->SetBranchStatus  ((m_branch_jet_truth + "_eta").c_str(), 1);
      tree->SetBranchAddress ((m_branch_jet_truth + "_eta").c_str(), &(ED_jet_truth->eta));

      tree->SetBranchStatus  ((m_branch_jet_truth + "_phi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_jet_truth + "_phi").c_str(), &(ED_jet_truth->phi));

      tree->SetBranchStatus  ((m_branch_jet_truth + "_E").c_str(), 1);
      tree->SetBranchAddress ((m_branch_jet_truth + "_E").c_str(), &(ED_jet_truth->E));
  }

  // L1 jets
  if ((CS->doTurnOns)&&(CS->useEmulation)){

      if(m_debug) Info("InputHandler()", "Setting L1 branches");

      tree->SetBranchStatus  ((m_branch_LVL1JetROI + "s_eta").c_str(), 1);
      tree->SetBranchAddress ((m_branch_LVL1JetROI + "s_eta").c_str(), &(L1D->eta));

      tree->SetBranchStatus  ((m_branch_LVL1JetROI + "s_phi").c_str(), 1);
      tree->SetBranchAddress ((m_branch_LVL1JetROI + "s_phi").c_str(), &(L1D->phi));

      tree->SetBranchStatus  ((m_branch_LVL1JetROI + "s_et8x8").c_str(), 1);
      tree->SetBranchAddress ((m_branch_LVL1JetROI + "s_et8x8").c_str(), &(L1D->Et));

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
    if(eventCounter == CS->doOnlyThisNumberOfEvents) return EL::StatusCode::SUCCESS;
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
  if ((CS->doOfflineJetKinematics)||(CS->doTriggerTruthResponse)||(CS->doTriggerOfflineResponse)||(CS->doTurnOns)){
      if(ED_jet->eta->size() < 2) return EL::StatusCode::SUCCESS;
  }
  // check trigger jets
  if ((CS->doTriggerJetKinematics)||(CS->doTriggerTruthResponse)||(CS->doTriggerOfflineResponse)||((CS->doTurnOns)&&(CS->useEmulation))){
      if(ED_trigJet->eta->size() < 2) return EL::StatusCode::SUCCESS;
  }
  // check truth jets
  if ((CS->doTruthJetKinematics)||(CS->doOfflineTruthResponse)||(CS->doTriggerTruthResponse)){
      if(ED_truthJet->eta->size() < 2) return EL::StatusCode::SUCCESS;
  }

  if (m_debugInExecute){
      cout << "eventCounter: " << eventCounter << endl;
  }

  // CHECK LArERROR ::: TODO: Change this into CUTS
  //  if (CS->checkLArError){
  // If we have LArERROR, skip this event
  //      if (m_LArError) return EL::StatusCode::SUCCESS;
  //  }

  /* ======= CUTS ======= */

  // Reset Cut Class
  cutH->Reset(TD->config_passedTriggers);
  cutH->Reset(TD->probe_passedTrigger);
  cutH->Reset(TD->ref_passedTrigger);

  // CHECK LArERROR
  if (CS->checkLArError) cutH->AddCut(!m_LArError); // if argument is false, event does NOT pass

  if (CS->doCuts){
      // Add simple cuts on observables
      //possible options:
      //"max": cutValue is maximal value
      //"min": cutValue is minimal value
      //"absmax": cutValue is maximal value of |observable|
      //"absmin": cutValue is minimal value of |observable|
      // cuts in GeV!

      //NOTE: Cuts are made on trigger jets!

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

  //TODO crosscheck
  /*  cout << "============================" << endl;
  cout << "The triggers that passes:..." << TD->event_passedTriggers->size() << endl;
  for (int n=0; n < TD->event_passedTriggers->size(); n++){
      cout << TD->event_passedTriggers->at(n) << endl;
      }*/

  /* ======= Matching ======= */

  // 1.1. Matching of trigger and offline jets

  // Initialise matchingIndex
  std::vector<int> matchingIndexList;
  std::vector<double> DeltaRList;


  // TODO Set Bool in CS for matching!
  /*  for (unsigned int i = 0; i < ED_trigJet->pt->size(); i++){
      matchingIndexList.push_back(-1);
      DeltaRList.push_back(-1);
      }*/
  // TODO RIGHT IF CONDITION FOR FINDBESTMATCHING
  //anaHandler->findBestMatching(ED_trigJet->pt, ED_trigJet->eta, ED_trigJet->phi, ED_trigJet->E, ED_jet->pt, ED_jet->eta, ED_jet->phi, ED_jet->E , matchingIndexList, DeltaRList, DeltaRMax);

  /* ======= CHECK JET BY JET => isGood ======= */

  if (m_debugInExecute) cout << "before isGood.." << endl;

  // OFFLINE JETS
  if ((CS->doOfflineJetKinematics)||(CS->doOfflineTruthResponse)||(CS->doTriggerOfflineResponse)){
      for (unsigned int i = 0; i < ED_jet->pt->size(); i++) {
        // first set true
        ED_jet->isGood[i] = true;
        if (CS->doCleaning) {
	      if (ED_jet->passedCleaning->at(i) == 0) ED_jet->isGood[i] = false;
        }
      }
  }

  // TRIGGER JETS
  if ((CS->doTriggerJetKinematics)||(CS->doTriggerTruthResponse)||(CS->doTriggerOfflineResponse)){
      for (unsigned int i = 0; i < ED_trigJet->pt->size(); i++) {
        // first set true
        ED_trigJet->isGood[i] = true;
        if (CS->doCleaning) {
	      if (ED_trigJet->passedCleaning->at(i) == 0) ED_trigJet->isGood[i] = false;
        }
      }
  }

  // TRUTH JETS
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
      m_map_triggerhistos[TD->config_triggerName.at(n)]->FillAll(ED_jet, ED_trigJet, ED_truthJet, ED_jet_truth, matchingIndexList, m_weight);
  }

  /* ======= Fill Turnons ======= */
  if (CS->doTurnOns){
      trigEffic->FillAll(TD, ED_jet, ED_trigJet, ED_truthJet, L1D, m_weight, CS);
  }


  eventCounter++;
  TD->UpdateCounting();
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
  trigEffic->DivideEfficiencyPlots(TD, CS);

  cout << "================== Results ====================" << endl;

  std::cout << "Number Events that passed the cuts: " << eventCounter << std::endl;

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

