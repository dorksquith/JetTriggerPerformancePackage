/*--------------------------------------------------------------------
KinematicMatrix.cxx

KinematicMatrix is a special class inherenting from HistogramMatrix that
contains all kinematic histograms of a specific jet type (offline, trigger, truth etc.)

created by Edgar Kellermann (edgar.kellermann@cern.ch)
--------------------------------------------------------------------*/

#include <JetTriggerPerformancePackage/KinematicMatrix.h>
#include <iostream>
#include <assert.h>


KinematicMatrix::KinematicMatrix(std::string key, std::string directory, int nthJet, ConfigStatus* a_CS):
    m_key(key), HistogramMatrix(directory, a_CS->doLeadSubleadThirdJetKinematics, a_CS->doNthJetKinematics), m_nthJet(nthJet), CS(a_CS)
{
    if (m_debug) std::cout << "Starting constructor KinematicMatrix()..." << std::endl;

    //Declare ToolsJTPP
    t = new ToolsJTPP();
}


KinematicMatrix::~KinematicMatrix()
{
  if (m_debug) std::cout << "Starting destructor KinematicMatrix()..." << std::endl;
}


void KinematicMatrix::Book(std::string name1, std::string name2, int nBins, double xMin, double xMax , EL::Worker* wk)
{
    if (m_debug) std::cout << "Starting Book():"<< name1 <<" ..." << std::endl;
    m_obs.push_back(name1);
    HistogramMatrix::Book(m_key + "_" + name1, m_key + "_" + name2, nBins, xMin, xMax , wk);


    //Additional Kinematic plots in bins of eta
    if (CS->doKinematicsInBinsOfEta){

	std::string stringAdd;
	// 1. region:
	stringAdd = "eta_0.0_" + t->ftos(etaLim1,2); //"eta_0.0_0.8";
	HistogramMatrix::Book(m_key + "_" + stringAdd + "_" + name1, m_key + "_" + stringAdd + "_" + name2, nBins, xMin, xMax , wk);
	// 2. region:
	stringAdd = "eta_" + t->ftos(etaLim1,2) + "_" + t->ftos(etaLim2,2); // "eta_0.8_1.6";
	HistogramMatrix::Book(m_key + "_" + stringAdd + "_" + name1, m_key + "_" + stringAdd + "_" + name2, nBins, xMin, xMax , wk);
	// 3. region:
	stringAdd = "eta_" + t->ftos(etaLim2,2) + "_" + t->ftos(etaLim3,2); // "eta_1.6_4.9";
	HistogramMatrix::Book(m_key + "_" + stringAdd + "_" + name1, m_key + "_" + stringAdd + "_" + name2, nBins, xMin, xMax , wk);
	// 4. region:
	stringAdd = "eta_" + t->ftos(etaLim3,2) + "_5.0"; // "eta_1.6_4.9";
	HistogramMatrix::Book(m_key + "_" + stringAdd + "_" + name1, m_key + "_" + stringAdd + "_" + name2, nBins, xMin, xMax , wk);
    }
}

void KinematicMatrix::Book(std::string name1, std::string name2, int nBins, double* xBinning, EL::Worker* wk)
{
    if (m_debug) std::cout << "Starting Book():"<< name1 <<" ..." << std::endl;
    m_obs.push_back(name1);
    HistogramMatrix::Book(m_key + "_" + name1, m_key + "_" + name2, nBins, xBinning, wk);


    //Additional Kinematic plots in bins of eta
    if (CS->doKinematicsInBinsOfEta){

	std::string stringAdd;
	// 1. region:
	stringAdd = "eta_0.0_" + t->ftos(etaLim1,2); //"eta_0.0_0.8";
	HistogramMatrix::Book(m_key + "_" + stringAdd + "_" + name1, m_key + "_" + stringAdd + "_" + name2, nBins, xBinning , wk);
	// 2. region:
	stringAdd = "eta_" + t->ftos(etaLim1,2) + "_" + t->ftos(etaLim2,2); // "eta_0.8_1.6";
	HistogramMatrix::Book(m_key + "_" + stringAdd + "_" + name1, m_key + "_" + stringAdd + "_" + name2, nBins, xBinning , wk);
	// 3. region:
	stringAdd = "eta_" + t->ftos(etaLim2,2) + "_" + t->ftos(etaLim3,2); // "eta_1.6_4.9";
	HistogramMatrix::Book(m_key + "_" + stringAdd + "_" + name1, m_key + "_" + stringAdd + "_" + name2,  nBins, xBinning , wk);
	// 4. region:
	stringAdd = "eta_" + t->ftos(etaLim3,2) + "_5.0"; // "eta_1.6_4.9";
	HistogramMatrix::Book(m_key + "_" + stringAdd + "_" + name1, m_key + "_" + stringAdd + "_" + name2,  nBins, xBinning , wk);
    }
}

void  KinematicMatrix::Fill(std::string obs, double value, double weight, bool isGood)
{
    if (m_debug) std::cout << "Starting Fill(): "<< obs <<" ..." << std::endl;
    HistogramMatrix::Fill(m_key + "_" + obs, value, isGood, weight);

}


void  KinematicMatrix::Fill(std::string obs, double value, std::vector<float>* etaVector, double weight, bool isGood)
{
    if (m_debug) std::cout << "Starting Fill(): "<< obs <<" ..." << std::endl;
    HistogramMatrix::Fill(m_key + "_" + obs, value, isGood, weight);

    //Additional Kinematic plots in bins of eta and check for cleaning
    if (CS->doKinematicsInBinsOfEta) {

	std::string stringAdd;
	// 1. region: 0 - 0.8
	if ((fabs(etaVector->at(0)) < etaLim1) && (fabs(etaVector->at(1)) < etaLim1)){
	    stringAdd = "eta_0.0_" + t->ftos(etaLim1,2);
	    HistogramMatrix::Fill(m_key + "_" + stringAdd + "_" + obs, value, isGood, weight);
	}
	// 2. region: 0.8 - 1.6
	if ((fabs(etaVector->at(0)) > etaLim1) && (fabs(etaVector->at(1)) > etaLim1) && (fabs(etaVector->at(0)) < etaLim2) && (fabs(etaVector->at(1)) < etaLim2)){
	    stringAdd = "eta_" + t->ftos(etaLim1,2) + "_" + t->ftos(etaLim2,2);
	    HistogramMatrix::Fill(m_key + "_" + stringAdd + "_" + obs, value, isGood, weight);
	}
	// 3. region: 1.6 - 4.9
	if ((fabs(etaVector->at(0)) > etaLim2) && (fabs(etaVector->at(1)) > etaLim2) && (fabs(etaVector->at(0)) < etaLim3) && (fabs(etaVector->at(1)) < etaLim3)){
	    stringAdd = "eta_" + t->ftos(etaLim2,2) + "_" + t->ftos(etaLim3,2);
	    HistogramMatrix::Fill(m_key + "_" + stringAdd + "_" + obs, value, isGood, weight);
	}
	// 3. region: 1.6 - 4.9
	if ((fabs(etaVector->at(0)) > etaLim3) && (fabs(etaVector->at(1)) > etaLim3) && (fabs(etaVector->at(0)) < 5.0) && (fabs(etaVector->at(1)) < 5.0)){
	    stringAdd = "eta_" + t->ftos(etaLim3,2) + "_5.0";
	    HistogramMatrix::Fill(m_key + "_" + stringAdd + "_" + obs, value, isGood, weight);
	}

    }

}


void  KinematicMatrix::Fill(std::string obs, int value, double weight, bool isGood)
{
    if (m_debug) std::cout << "Starting Fill(): "<< obs <<" ..." << std::endl;
    HistogramMatrix::Fill(m_key + "_" + obs, value, isGood, weight);

}


void  KinematicMatrix::Fill(std::string obs, std::vector<float>* valueVector, double weight, std::vector<bool>* isGood)
{
    if (m_debug) std::cout << "Starting Fill(): "<< obs <<" ..." << std::endl;
    HistogramMatrix::Fill(m_key + "_" + obs, valueVector, isGood, weight);

}


void  KinematicMatrix::Fill(std::string obs, std::vector<float>* valueVector, std::vector<float>* etaVector, double weight, std::vector<bool>* isGood)
{
    if (m_debug) std::cout << "Starting Fill(): "<< obs <<" ..." << std::endl;
    HistogramMatrix::Fill(m_key + "_" + obs, valueVector, isGood, weight);


    //Additional Kinematic plots in bins of eta
    if (CS->doKinematicsInBinsOfEta){ // A MESS OF A CODE!
	if (m_debug) std::cout << "doKinematicsInBinsOfEta..." << std::endl;

	for (unsigned int n=0; n < valueVector->size(); n++){

	    std::string stringAdd;
	    // 1. region: 0 - 0.8
	    if (fabs(etaVector->at(n)) < etaLim1){
		stringAdd = "eta_0.0_" + t->ftos(etaLim1,2);
		if (m_debug) std::cout << "doKinematicsInBinsOfEta..." << stringAdd << std::endl;
		HistogramMatrix::FillConsiderLeadSubThird(m_key + "_" + stringAdd + "_" + obs, valueVector->at(n), n, isGood->at(n), weight);
	    }
	    // 2. region: 0.8 - 1.6
	    if ((fabs(etaVector->at(n)) > etaLim1) && (fabs(etaVector->at(n)) < etaLim2)){
		stringAdd = "eta_" + t->ftos(etaLim1,2) + "_" + t->ftos(etaLim2,2);
		if (m_debug) std::cout << "doKinematicsInBinsOfEta..." << stringAdd << std::endl;
		HistogramMatrix::FillConsiderLeadSubThird(m_key + "_" + stringAdd + "_" + obs, valueVector->at(n), n, isGood->at(n), weight);
	    }
	    // 3. region: 1.6 - 4.9
	    if ((fabs(etaVector->at(n)) > etaLim2) && (fabs(etaVector->at(n)) < etaLim3)){
		stringAdd = "eta_" + t->ftos(etaLim2,2) + "_" + t->ftos(etaLim3,2);
		if (m_debug) std::cout << "doKinematicsInBinsOfEta..." << stringAdd << std::endl;
		HistogramMatrix::FillConsiderLeadSubThird(m_key + "_" + stringAdd + "_" + obs, valueVector->at(n), n, isGood->at(n), weight);
	    }
	    // 4. region: 4.9 - 5.0
	    if ((fabs(etaVector->at(n)) > etaLim3) && (fabs(etaVector->at(n)) < 5.0)){
		stringAdd = "eta_" + t->ftos(etaLim3,2) + "_5.0";
		if (m_debug) std::cout << "doKinematicsInBinsOfEta..." << stringAdd << std::endl;
		HistogramMatrix::FillConsiderLeadSubThird(m_key + "_" + stringAdd + "_" + obs, valueVector->at(n), n, isGood->at(n), weight);
	    }

	}
    }

}


void KinematicMatrix::FillAll(EventData* ED, double weight)
{
   if (m_debug) std::cout << "Starting FillAll()..." << std::endl;

   this->Fill("E", ED->E, ED->eta, weight, &(ED->isGood));
   this->Fill("pt", ED->pt, ED->eta, weight, &(ED->isGood));
   this->Fill("phi", ED->phi, ED->eta, weight, &(ED->isGood));
   this->Fill("eta", ED->eta, ED->eta, weight, &(ED->isGood));

   //more observables
   if(CS->doMjj)         this->Fill("mjj", ED->Mjj(), ED->eta, weight, ED->isObsGood(0,1));
   if(CS->doM23)         this->Fill("m23", ED->M23(), ED->eta, weight, ED->isObsGood(1,2));
   if(CS->doyStar)       this->Fill("yStar", ED->YStar(), ED->eta, weight, ED->isObsGood(0,1));
   if(CS->doDeltaPhi)    this->Fill("deltaPhi", ED->deltaPhi, ED->eta, weight, ED->isObsGood(0,1));
   if(CS->doPTBalance)   this->Fill("pTBalance", ED->pTBalance, ED->eta, weight, true); // NOTE: we dont consider cleaning for the following observ.
   if(CS->doHT)          this->Fill("HT", ED->GetHT(), ED->eta, weight, true);
   if(CS->doMHT)         this->Fill("MHT", ED->MHT, ED->eta, weight, true);
   if(CS->doMHTPhi)      this->Fill("MHTPhi", ED->MHTPhi, ED->eta, weight, true);
   if(CS->doEMFrac)      this->Fill("EMFrac", ED->EMFrac, ED->eta, weight, &(ED->isGood));
   if(CS->doHECFrac)     this->Fill("HECFrac", ED->HECFrac, ED->eta, weight, &(ED->isGood));
   if(CS->doFracSamplingMax) this->Fill("FracSamplingMax", ED->FracSamplingMax, ED->eta, weight, &(ED->isGood));

}


