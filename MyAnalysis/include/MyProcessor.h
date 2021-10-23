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

class MyProcessor : public Processor
{
  public:
    virtual MyProcessor* newProcessor() { return new MyProcessor; }

    MyProcessor();

    virtual void init();
    virtual void processEvent(LCEvent* evt);
    virtual void end();

    MyProcessor(const MyProcessor& toCopy) = delete;
    void operator=(const MyProcessor& toCopy) = delete;

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

    float totalEnergy = 0;
    
    float mc_GammaEnergy = 0;
    float mc_GammaPt = 0;
    float mc_GammaMass = 0;
    
    float collEnergy = 0;
    float collPt = 0;
    
    int mc_NuType = 0;
    float mc_NuEnergy = 0;
    float mc_NuPt = 0;
    float mc_NuMass = 0;

    float mc_HiggsEnergy = 0;
    float mc_HiggsPt = 0;
    float mc_HiggsMass = 0;
    
    int mc_HiggsDecay = 0;
    float mc_Decay1Energy = 0;
    float mc_Decay1Pt = 0;
    float mc_Decay1Mass = 0;
    float mc_Decay2Energy = 0;
    float mc_Decay2Pt = 0;
    float mc_Decay2Mass = 0;
    float mc_Cos = 0;
    
    int mc_DecayProcess = 0;
    int mc_DecayNuProcess = 0;
    
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
