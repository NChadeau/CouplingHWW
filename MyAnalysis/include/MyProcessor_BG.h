#pragma once

#include "lcio.h"
#include "marlin/Processor.h"

class TFile;
class TTree;
namespace IMPL
{
class ReconstructedParticleImpl;
class MCParticleImpl;
} // namespace IMPL

using namespace lcio;
using namespace marlin;

class MyProcessor_BG : public Processor
{
  public:
    virtual MyProcessor_BG* newProcessor() { return new MyProcessor_BG; }

    MyProcessor_BG();

    virtual void init();
    virtual void processEvent(LCEvent* evt);
    virtual void end();

    MyProcessor_BG(const MyProcessor_BG& toCopy) = delete;
    void operator=(const MyProcessor_BG& toCopy) = delete;

  private:
    void clear();

  private:
    int nEventsProcessed = 0;

    std::string mcParticleCollectionName = "";
    std::string reconstructedParticleCollectionName = "";

    std::vector<IMPL::ReconstructedParticleImpl*> reconstructedParticles{};

    std::string rootFileName = "";
    TFile*      file = nullptr;
    TTree*      tree = nullptr;
    
    double y_cut = 0.002;
    int number_jets = 4;
    
    int reco_njets = 0;
    
    float reco_Y12 = 0;
    float reco_Y13 = 0;
    float reco_Y14 = 0;
    float reco_Y23 = 0;
    float reco_Y24 = 0;
    float reco_Y34 = 0;
    
    float reco_E4jets = 0;
    float reco_Pt4jets = 0;
    float reco_M4jets = 0;
    
    float reco_WBigMass_Energy = 0;
    float reco_WBigMass_Pt = 0;
    float reco_WBigMass_Mass = 0;
    float reco_WBigMass_CosJets = 0;
    float reco_WSmallMass_Energy = 0;
    float reco_WSmallMass_Pt = 0;
    float reco_WSmallMass_Mass = 0;
    float reco_WSmallMass_CosJets = 0;
    float reco_Cos = 0;
    
    int reco_BoolDumpEvent = 0;
    
};
