/*--------------------------------------------------------------------
EventData.cxx

EventData contains the event observables of one jet type (offline, trigger, truth etc.).

created by Edgar Kellermann (edgar.kellermann@cern.ch)
--------------------------------------------------------------------*/


#include <JetTriggerPerformancePackage/EventData.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include "TLorentzVector.h"



EventData::EventData(std::string jetType, ConfigStatus* a_CS):
    CS(a_CS), m_jetType(jetType), E(nullptr), pt(nullptr), eta(nullptr), phi(nullptr), mjj(0), m23(0), passedCleaning(nullptr), Timing(nullptr)
{
  if (m_debug) std::cout << "Starting constructor EventData()..." << std::endl;
}

EventData::~EventData()
{
  if (m_debug) std::cout << "Starting destructor EventData()..." << std::endl;
}

bool EventData::isObsGood(int i, int j)
{
    if (m_debug) std::cout << "Starting destructor isLSGood()..." << std::endl;
    if ((isGood.at(i)) && (isGood.at(j))) return true;
    else return false;
}

// calculate transverse Energy
float EventData::GetEt(int n)
{
    if (m_debug) std::cout << "Starting destructor GetEt()..." << std::endl;

    TLorentzVector vec;
    vec.SetPtEtaPhiE(pt->at(n), eta->at(n), phi->at(n), E->at(n));

    return vec.Et();
}

float EventData::GetMjj(int i=0, int j=1)
{
    if (m_debug) std::cout << "Starting method Mjj()..." << std::endl;

    // crosscheck vector size
    if ((pt->size() < (i+1)) || (pt->size() < (j+1))) return 0.0;

    TLorentzVector vec1;
    TLorentzVector vec2;

    // fill Lorentz vectors
    vec1.SetPtEtaPhiE(pt->at(i), eta->at(i), phi->at(i), E->at(i));
    vec2.SetPtEtaPhiE(pt->at(j), eta->at(j), phi->at(j), E->at(j));

    // combine vectors
    TLorentzVector vecSum = vec1 + vec2;

    // return invariant mass
    return vecSum.M();
}

float EventData::GetYStar(int i=0, int j=1)
{
    if (m_debug) std::cout << "Starting method yStar()..." << std::endl;

    // crosscheck vector size
    if ((pt->size() < (i+1)) || (pt->size() < (j+1))) return -999.0;

    TLorentzVector vec1;
    TLorentzVector vec2;

    // fill Lorentz vectors
    vec1.SetPtEtaPhiE(pt->at(i), eta->at(i), phi->at(i), E->at(i));
    vec2.SetPtEtaPhiE(pt->at(j), eta->at(j), phi->at(j), E->at(j));

    // calculate rapidity and return yStar
    return 0.5* fabs(vec1.Rapidity() - vec2.Rapidity());
}

float EventData::GetHT()
{
    if (m_debug) std::cout << "Starting method GetHT()..." << std::endl;

    float HT = 0.0;

    // sum over all pt (if pt > 50 GeV and |eta| < 2.8)
    for (int j=0; j < pt->size(); j++){
        if ((pt->at(j) > 50.0) && (abs(eta->at(j)) < 2.8)) HT += pt->at(j);
    }

    return HT;
}

float EventData::GetHT_triggerLevel()
{
    if (m_debug) std::cout << "Starting method GetHT_triggerLevel()..." << std::endl;

    float HT = 0.0;
    float Et;

    // sum over all pt (if pt > 30 GeV and |eta| < 3.2)
    for (int j=0; j < pt->size(); j++){
        Et = this->GetEt(j);
        if (( Et > 30.0) && (fabs(eta->at(j)) < 3.2)) HT += Et;
    }

    return HT;
}

float EventData::Mjj()
{
    if (m_debug) std::cout << "Starting method Mjj()..." << std::endl;

    if (CS->calculateMjj) return this->GetMjj(0,1);
    else return mjj;
}

float EventData::M23()
{
    if (m_debug) std::cout << "Starting method M23()..." << std::endl;

    if (CS->calculateMjj) return this->GetMjj(1,2);
    else return m23;
}

float EventData::YStar()
{
    if (m_debug) std::cout << "Starting method YStar()..." << std::endl;

    if (CS->calculateYStar) return this->GetYStar(0,1);
    else return yStar;
}
