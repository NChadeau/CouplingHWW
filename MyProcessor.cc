#include "MyProcessor.h"

#include <iostream>
#include <cmath>

#include <EVENT/LCEvent.h>
#include <TFile.h>
#include <TTree.h>

#include <EVENT/LCCollection.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <fastjet/ClusterSequence.hh>

#include <CLHEP/Vector/ThreeVector.h>

using namespace fastjet;
using namespace std;

MyProcessor aMyProcessor;

double jet_parameter(const fastjet::PseudoJet& jet_1, const fastjet::PseudoJet& jet_2, double total_E) {   
    CLHEP::Hep3Vector p1(jet_1.px(), jet_1.py(), jet_1.pz());
    CLHEP::Hep3Vector p2(jet_2.px(), jet_2.py(), jet_2.pz());
    double log_Y = (-1)*std::log(jet_1.E()*jet_2.E() * (1 - (p1.dot(p2))/(p1.mag()*p2.mag()))* 2./(total_E*total_E));  
    return log_Y;
}

double dist_mass(const fastjet::PseudoJet& jet_1, const fastjet::PseudoJet& jet_2) {
    double d = std::min(std::abs(jet_1.m() - 80.379) + std::abs(jet_2.m() - 37.654), std::abs(jet_2.m() - 80.379) + std::abs(jet_1.m() - 37.654));
    return d;
}

bool test_order(const fastjet::PseudoJet& jet_1, const fastjet::PseudoJet& jet_2, const double d) {
    if (d == std::abs(jet_1.m() - 80.379) + std::abs(jet_2.m() - 37.654)) {
        bool order = true;
        return order;
    } else {
        bool order = false;
        return order;
    }
}

double cos(const fastjet::PseudoJet& jet_1, const fastjet::PseudoJet& jet_2) {
    double cos = (jet_1.px()*jet_2.px() + jet_1.py()*jet_2.py() + jet_1.pz()*jet_2.pz())/(jet_1.modp()*jet_2.modp());
    return cos;
}

MyProcessor::MyProcessor()
    : Processor("MyProcessor")
{
    registerProcessorParameter("RootFileName", "File name for the root output", rootFileName, std::string("test.root"));

    registerProcessorParameter("MCParticlesCollectionName", "Name of the MC particles collection", mcParticleCollectionName, std::string("MCParticlesSkimmed"));
    
    registerProcessorParameter("ReconstructedParticlesCollectionName", "Name of the reconstructed particles collection",reconstructedParticleCollectionName, std::string("PandoraPFOs"));
}

void MyProcessor::init() {
    std::cout << "MyProcessor::init()" << std::endl;
    file = new TFile(rootFileName.c_str(), "RECREATE");
    tree = new TTree("tree", "tree");

    //Total energy
    tree->Branch("totalEnergy", &totalEnergy);
    
    //Energy and Pt of radiations from acceleration
    tree->Branch("mc_GammaEnergy", &mc_GammaEnergy);
    tree->Branch("mc_GammaPt", &mc_GammaPt);
    tree->Branch("mc_GammaMass", &mc_GammaMass);
    
    //Energy and Pt of the collision
    tree->Branch("collEnergy", &collEnergy);
    tree->Branch("collPt", &collPt);
    
    //Energy, Pt and mass of the Higgs
    tree->Branch("mc_HiggsEnergy", &mc_HiggsEnergy);
    tree->Branch("mc_HiggsPt" , &mc_HiggsPt);
    tree->Branch("mc_HiggsMass", &mc_HiggsMass);
    
    //Type, energy, Pt and invariant mass of the neutrinos
    tree->Branch("mc_NuType", &mc_NuType);
    tree->Branch("mc_NuEnergy", &mc_NuEnergy);
    tree->Branch("mc_NuPt", &mc_NuPt);
    tree->Branch("mc_NuMass", &mc_NuMass);
    
    //Type, Energy, Pt and mass of the each Higgs' decays and the angle between them
    tree->Branch("mc_HiggsDecay", &mc_HiggsDecay);
    tree->Branch("mc_Decay1Energy", &mc_Decay1Energy);
    tree->Branch("mc_Decay1Pt" , &mc_Decay1Pt);
    tree->Branch("mc_Decay1Mass", &mc_Decay1Mass);
    tree->Branch("mc_Decay2Energy", &mc_Decay2Energy);
    tree->Branch("mc_Decay2Pt" , &mc_Decay2Pt);
    tree->Branch("mc_Decay2Mass", &mc_Decay2Mass);
    tree->Branch("mc_Cos" , &mc_Cos);
    
    //Type of decay: hadronic, leptinic or semileptonic and the number of neutrini
    tree->Branch("mc_DecayProcess", &mc_DecayProcess);
    tree->Branch("mc_DecayNuProcess", &mc_DecayNuProcess);
    
    //Number of jets
    tree->Branch("reco_njets", &reco_njets);
     
    //Reconstructed jets parameters
    tree->Branch("reco_Y12", &reco_Y12);
    tree->Branch("reco_Y13", &reco_Y13);
    tree->Branch("reco_Y14", &reco_Y14);
    tree->Branch("reco_Y23", &reco_Y23);
    tree->Branch("reco_Y24", &reco_Y24);
    tree->Branch("reco_Y34", &reco_Y34);
    
    //Energy, impulsion and invariant mass of the 4 jets
    tree->Branch("reco_E4jets", &reco_E4jets);
    tree->Branch("reco_Pt4jets", &reco_Pt4jets);
    tree->Branch("reco_M4jets", &reco_M4jets);
    
    //Energy, impulsion and mass of bi-jets corresponding with mc_HiggsDecay == 2424 (M == 80.379, m == 37.654)
    tree->Branch("reco_WBigMass_Energy", &reco_WBigMass_Energy);
    tree->Branch("reco_WBigMass_Pt", &reco_WBigMass_Pt);
    tree->Branch("reco_WBigMass_Mass", &reco_WBigMass_Mass);
    tree->Branch("reco_WBigMass_CosJets", &reco_WBigMass_CosJets);
    tree->Branch("reco_WSmallMass_Energy", &reco_WSmallMass_Energy);
    tree->Branch("reco_WSmallMass_Pt", &reco_WSmallMass_Pt);
    tree->Branch("reco_WSmallMass_Mass", &reco_WSmallMass_Mass);
    tree->Branch("reco_WSmallMass_CosJets", &reco_WSmallMass_CosJets);
    tree->Branch("reco_Cos", &reco_Cos);
    
    //To know if there is at least 4 reconstructed particles or not
    tree->Branch("reco_BoolDumpEvent", &reco_BoolDumpEvent);
    
}

void MyProcessor::clear() {
    reconstructedParticles.clear();
    totalEnergy = 0;
}

void MyProcessor::processEvent(LCEvent* evt) {
    clear();

    LCCollection* reconstructedParticleCollection = evt->getCollection(reconstructedParticleCollectionName);
    const int     nRecoParticles = reconstructedParticleCollection->getNumberOfElements();
    for (int index = 0; index < nRecoParticles; ++index) {
        auto particle = dynamic_cast<IMPL::ReconstructedParticleImpl*>(reconstructedParticleCollection->getElementAt(index));
        reconstructedParticles.push_back(particle);
    }

    LCCollection* mcParticleCollection = evt->getCollection(mcParticleCollectionName);
    
    auto gamma0 = dynamic_cast<IMPL::MCParticleImpl*>(mcParticleCollection->getElementAt(6));
    auto gamma1 = dynamic_cast<IMPL::MCParticleImpl*>(mcParticleCollection->getElementAt(7));
    mc_GammaEnergy = gamma0->getEnergy() + gamma1->getEnergy();
    const double gamma_px = gamma0->getMomentum()[0] + gamma1->getMomentum()[0];
    const double gamma_py = gamma0->getMomentum()[1] + gamma1->getMomentum()[1];
    const double gamma_pz = gamma0->getMomentum()[2] + gamma1->getMomentum()[2];
    mc_GammaPt = std::sqrt(gamma_px * gamma_px + gamma_py * gamma_py);
    double test = mc_GammaEnergy*mc_GammaEnergy - (mc_GammaPt*mc_GammaPt + gamma_pz*gamma_pz);
    if (test < 0) {
        mc_GammaMass = -(1)*std::sqrt((-1)*test);
    } else {
        mc_GammaMass = std::sqrt(test);
    }
    
    auto nu0 = dynamic_cast<IMPL::MCParticleImpl*>(mcParticleCollection->getElementAt(8));
    auto nu1 = dynamic_cast<IMPL::MCParticleImpl*>(mcParticleCollection->getElementAt(9));
    mc_NuType = std::abs(nu0->getPDG());
    mc_NuEnergy = nu0->getEnergy() + nu1->getEnergy();
    const double nu_px = nu0->getMomentum()[0] + nu1->getMomentum()[0];
    const double nu_py = nu0->getMomentum()[1] + nu1->getMomentum()[1];
    const double nu_pz = nu0->getMomentum()[2] + nu1->getMomentum()[2];
    mc_NuPt = std::sqrt(nu_px * nu_px + nu_py * nu_py);
    mc_NuMass = std::sqrt(mc_NuEnergy * mc_NuEnergy - mc_NuPt * mc_NuPt - nu_pz * nu_pz);

    auto higgsMCParticle = dynamic_cast<IMPL::MCParticleImpl*>(mcParticleCollection->getElementAt(10));
    mc_HiggsEnergy = higgsMCParticle->getEnergy();
    const double higgs_px = higgsMCParticle->getMomentum()[0];
    const double higgs_py = higgsMCParticle->getMomentum()[1];
    const double higgs_pz = higgsMCParticle->getMomentum()[2];
    mc_HiggsPt = std::sqrt(higgs_px * higgs_px + higgs_py * higgs_py);
    mc_HiggsMass = std::sqrt(mc_HiggsEnergy * mc_HiggsEnergy - mc_HiggsPt * mc_HiggsPt - higgs_pz * higgs_pz);
    
    totalEnergy = mc_HiggsEnergy + mc_GammaEnergy + mc_NuEnergy;
    
    collEnergy = mc_HiggsEnergy + mc_NuEnergy;
    collPt = (nu_px + higgs_px) * (nu_px + higgs_px) + (nu_py + higgs_py) * (nu_px + higgs_px);
    
    auto decay1 = dynamic_cast<IMPL::MCParticleImpl*>(mcParticleCollection->getElementAt(11));
    auto decay2 = dynamic_cast<IMPL::MCParticleImpl*>(mcParticleCollection->getElementAt(12));
    // Higgs into 2 gluons
    if (std::abs(decay1->getPDG()) == 21 && std::abs(decay2->getPDG()) == 21) {
        mc_HiggsDecay = 2121; //2 gluons
        mc_DecayProcess = 808; //hadronic
        mc_DecayNuProcess = 0; //0 neutrino
    // Higgs into 2 photons
    } else if (std::abs(decay1->getPDG()) == 22 && std::abs(decay2->getPDG()) == 22) {
        mc_HiggsDecay = 2222; //2 photons
        mc_DecayProcess = 1616; //2 photons
        mc_DecayNuProcess = 0; //0 neutrino
    //Higgs into 1 photon and 1 Z boson
    } else if ((std::abs(decay1->getPDG()) == 22 && std::abs(decay2->getPDG()) == 23) || (std::abs(decay1->getPDG()) == 23 && std::abs(decay2->getPDG()) == 22)) {
        mc_HiggsDecay = 2223; //1 photon and 1 Z boson
        if (std::abs(decay1->getPDG()) == 23) {
            double d1 = std::abs(decay1->getDaughters()[0]->getPDG());
            if (d1 < 7) {
            	mc_DecayProcess = 1608; //photon and hadronic
            	mc_DecayNuProcess = 0; //0 neutrino
            } else {
            	mc_DecayProcess = 1612; //photon and leptonic
            	if (d1 == 12 || d1 == 14 || d1 == 16) {
            	    mc_DecayNuProcess = 2; //2 neutrinos
            	} else {
            	    mc_DecayNuProcess = 0; //0 neutrino
            	}
            }
        } else {
            double d2 = std::abs(decay2->getDaughters()[0]->getPDG());
            if (d2 < 7) {
            	mc_DecayProcess = 1608; //photon and hadronic
            	mc_DecayNuProcess = 0; //0 neutrino
            } else {
            	mc_DecayProcess = 1612; //photon and leptonic
            	if (d2 == 12 || d2 == 14 || d2 == 16) {
            	    mc_DecayNuProcess = 2; //2 neutrinos
            	} else {
            	    mc_DecayNuProcess = 0; //0 neutrino
            	}
            }
        }
    //Higgs into 2 Z bosons
    } else if (std::abs(decay1->getPDG()) == 23 && std::abs(decay2->getPDG()) == 23) {
        mc_HiggsDecay = 2323; //2 Z bosons
        double d1 = std::abs(decay1->getDaughters()[0]->getPDG());
        double d2 = std::abs(decay2->getDaughters()[0]->getPDG());
        if (d1 < 7 && d2 < 7) {
            mc_DecayProcess = 808; //hadronic
            mc_DecayNuProcess = 0; //0 neutrino
        } else if ((d1 < 7 && d2 > 10 && d2 < 17) || (d2 < 7 && d1 > 10 && d1 < 17)) {
            mc_DecayProcess = 1208; //semileptonic
            if (d1 < 7 && (d2 == 12 || d2 == 14 || d2 == 16)) {
            	mc_DecayNuProcess = 2; //2 neutrinos
            } else if (d2 < 7 && (d1 == 12 || d1 == 14 || d1 == 16)) {
                mc_DecayNuProcess = 2; //2 neutrinos
            } else {
                mc_DecayNuProcess = 0; //0 neutrino
            }
        } else {
            mc_DecayProcess = 1212; //leptonic 
            if ((d1 == 12 || d1 == 14 || d1 == 16) && (d2 == 12 || d2 == 14 || d2 == 16)) {
            	mc_DecayNuProcess = 4; //4 neutrinos
            } else if ((d1 == 12 || d1 == 14 || d1 == 16) || (d2 == 12 || d2 == 14 || d2 == 16)) {
                mc_DecayNuProcess = 2; //2 neutrinos
            } else {
                mc_DecayNuProcess = 0; //0 neutrino
            }
        }
    // Higgs into 2 W bosons
    } else if (std::abs(decay1->getPDG()) == 24 && std::abs(decay2->getPDG()) == 24) {
        mc_HiggsDecay = 2424; //2 W bosons
        double d1 = std::abs(decay1->getDaughters()[0]->getPDG());
        double d2 = std::abs(decay2->getDaughters()[0]->getPDG());
        if (d1 < 7 && d2 < 7) {
            mc_DecayProcess = 808; //hadronic
            mc_DecayNuProcess = 0; //0 neutrino
        } else if ((d1 < 7 && (d2 < 17 && d2 > 10)) || (d2 < 7 && (d1 < 17 && d1 > 10))) {
            mc_DecayProcess = 1208; //semileptonic
            mc_DecayNuProcess = 1; //1 neutrino
        } else {
            mc_DecayProcess = 1212; //leptonic
            mc_DecayNuProcess = 2; //2 neutrinos
        }
    //Higgs into 2 fermions
    } else {
    	double d1 = std::abs(decay1->getDaughters()[0]->getPDG());
        mc_HiggsDecay = d1*10 + d1; 
        if (d1 < 7) {
            mc_DecayProcess = 808; //hadronic
            mc_DecayNuProcess = 0; //0 neutrino
        } else {
            mc_DecayProcess = 1212; //leptonic
            if (d1 == 12 || d1 == 14 || d1 == 16) {
                mc_DecayNuProcess = 2; //2 neutrinos
            } else {
                mc_DecayNuProcess = 0; //0 neutrino
            }
        }
    }
    double e1 = decay1->getEnergy();
    double px_1 = decay1->getMomentum()[0]; 
    double py_1 = decay1->getMomentum()[1];
    double pz_1 = decay1->getMomentum()[2];
    double pt_1 = std::sqrt(px_1 * px_1 + py_1 * py_1);
    double p_1 = std::sqrt(px_1 * px_1 + py_1 * py_1 + pz_1 * pz_1);
    double e2 = decay2->getEnergy();
    double px_2 = decay2->getMomentum()[0]; 
    double py_2 = decay2->getMomentum()[1];
    double pz_2 = decay2->getMomentum()[2];
    double pt_2 = std::sqrt(px_2 * px_2 + py_2 * py_2);
    double p_2 = std::sqrt(px_2 * px_2 + py_2 * py_2 + pz_2 * pz_2);
    //Sort them such that mc_Decay1Energy > mc_Decay2Energy
    if (e1 > e2) {
        mc_Decay1Energy = e1;
        mc_Decay1Pt = pt_1;
        mc_Decay1Mass = std::sqrt(e1 * e1 - p_1 * p_1);
        mc_Decay2Energy = e2;
        mc_Decay2Pt = pt_2;
        mc_Decay2Mass = std::sqrt(e2 * e2 - p_2 * p_2);
    } else {
        mc_Decay1Energy = e2;
        mc_Decay1Pt = pt_2;
        mc_Decay1Mass = std::sqrt(e2 * e2 - p_2 * p_2);
        mc_Decay2Energy = e1;
        mc_Decay2Pt = pt_1;
        mc_Decay2Mass = std::sqrt(e1 * e1 - p_1 * p_1);
    }
    mc_Cos = (px_1 * px_2 + py_1 * py_2 + pz_1 * pz_2) / (p_1 * p_2);
    
    //Jets study
    std::vector<fastjet::PseudoJet> Particles ;
    for (int index = 0; index < nRecoParticles; ++index) {
        auto particle = dynamic_cast<IMPL::ReconstructedParticleImpl*>(reconstructedParticleCollection->getElementAt(index));
        fastjet::PseudoJet jets_init(particle->getMomentum()[0], particle->getMomentum()[1], particle->getMomentum()[2], particle->getEnergy());
        Particles.push_back(jets_init);
    }
    fastjet::JetDefinition j(fastjet::ee_kt_algorithm);
    fastjet::ClusterSequence cs(Particles, j);
    std::vector<fastjet::PseudoJet> jets_final = cs.exclusive_jets_ycut(y_cut);
    reco_njets = jets_final.size();
    
    if (nRecoParticles > 3) { // we need at least 4 reconstructed particles to study Y_ij distributions
        reco_BoolDumpEvent = 0;
        std::vector<fastjet::PseudoJet> Particles_bis ;
        for (int index = 0; index < nRecoParticles; ++index) {
            auto particle = dynamic_cast<IMPL::ReconstructedParticleImpl*>(reconstructedParticleCollection->getElementAt(index));
            fastjet::PseudoJet jets_init(particle->getMomentum()[0], particle->getMomentum()[1], particle->getMomentum()[2], particle->getEnergy());
            Particles_bis.push_back(jets_init);
        }
        fastjet::JetDefinition j_bis(fastjet::ee_kt_algorithm);
        fastjet::ClusterSequence cs_bis(Particles_bis, j_bis);
        std::vector<fastjet::PseudoJet> jets_final_bis = sorted_by_E(cs_bis.exclusive_jets(number_jets));
        fastjet::PseudoJet jets_total = jets_final_bis[0] + jets_final_bis[1] + jets_final_bis[2] + jets_final_bis[3];
        reco_Y12 = jet_parameter(jets_final_bis[0], jets_final_bis[1], 250);
        reco_Y13 = jet_parameter(jets_final_bis[0], jets_final_bis[2], 250);
        reco_Y14 = jet_parameter(jets_final_bis[0], jets_final_bis[3], 250);
        reco_Y23 = jet_parameter(jets_final_bis[1], jets_final_bis[2], 250);
        reco_Y24 = jet_parameter(jets_final_bis[1], jets_final_bis[3], 250);
        reco_Y34 = jet_parameter(jets_final_bis[2], jets_final_bis[3], 250);
        reco_E4jets = jets_total.E();
        reco_Pt4jets = jets_total.pt();
        reco_M4jets = jets_total.m();
        std::vector<std::vector<fastjet::PseudoJet>> jets;
        std::vector<fastjet::PseudoJet> jets_1 = {jets_final_bis[0] + jets_final_bis[1], jets_final_bis[2] + jets_final_bis[3]};
        std::vector<fastjet::PseudoJet> jets_2 = {jets_final_bis[0] + jets_final_bis[2], jets_final_bis[1] + jets_final_bis[3]};
        std::vector<fastjet::PseudoJet> jets_3 = {jets_final_bis[0] + jets_final_bis[3], jets_final_bis[1] + jets_final_bis[2]};
        jets.push_back(jets_1);
        jets.push_back(jets_2);
        jets.push_back(jets_3);
        std::vector<std::vector<double>> cos_jets;
        std::vector<double> cos_jets_1 = {cos(jets_final_bis[0], jets_final_bis[1]), cos(jets_final_bis[2], jets_final_bis[3])};
        std::vector<double> cos_jets_2 = {cos(jets_final_bis[0], jets_final_bis[2]), cos(jets_final_bis[1], jets_final_bis[3])};
        std::vector<double> cos_jets_3 = {cos(jets_final_bis[0], jets_final_bis[3]), cos(jets_final_bis[1], jets_final_bis[2])};
        cos_jets.push_back(cos_jets_1);
        cos_jets.push_back(cos_jets_2);
        cos_jets.push_back(cos_jets_3);
        double d = 1000;
        int order = false;
        fastjet::PseudoJet jet_init_1 = {0, 0, 0, 0};
        fastjet::PseudoJet jet_init_2 = {0, 0, 0, 0};
        std::vector<fastjet::PseudoJet> selected_jets = {jet_init_1, jet_init_2};
        std::vector<double> selected_cos_jets = {0, 0};
        for (int i = 0; i < 3; i++) {
            if (dist_mass(jets[i][0], jets[i][2]) < d) {
                d = dist_mass(jets[i][0], jets[i][2]);
                order = test_order(jets[i][0], jets[i][2], d);
                selected_jets = jets[i];
                selected_cos_jets = cos_jets[i];
            }
        } 
        if (order == true) {
            reco_WBigMass_Energy = selected_jets[0].E();
            reco_WBigMass_Pt = selected_jets[0].pt();
            reco_WBigMass_Mass = selected_jets[0].m();
            reco_WBigMass_CosJets = selected_cos_jets[0];
            reco_WSmallMass_Energy = selected_jets[1].E();
            reco_WSmallMass_Pt = selected_jets[1].pt();
            reco_WSmallMass_Mass = selected_jets[1].m();
            reco_WSmallMass_CosJets = selected_cos_jets[1];
        } else {
           reco_WBigMass_Energy = selected_jets[1].E();
           reco_WBigMass_Pt = selected_jets[1].pt();
           reco_WBigMass_Mass = selected_jets[1].m();
           reco_WBigMass_CosJets = selected_cos_jets[1];
           reco_WSmallMass_Energy = selected_jets[0].E();
           reco_WSmallMass_Pt = selected_jets[0].pt();
           reco_WSmallMass_Mass = selected_jets[0].m();
           reco_WSmallMass_CosJets = selected_cos_jets[0];
        }
        reco_Cos = cos(selected_jets[0], selected_jets[1]);
    } else {
        reco_BoolDumpEvent = 1;
        reco_Y12 = 1000;
        reco_Y13 = 1000;
        reco_Y14 = 1000;
        reco_Y23 = 1000;
        reco_Y24 = 1000;
        reco_Y34 = 1000;
        reco_E4jets = 1000;
        reco_Pt4jets = 1000;
        reco_M4jets = 1000;
        reco_WBigMass_Energy = 1000;
        reco_WBigMass_Pt = 1000;
        reco_WBigMass_Mass = 1000;
        reco_WBigMass_CosJets = 1000;
        reco_WSmallMass_Energy = 1000;
        reco_WSmallMass_Pt = 1000;
        reco_WSmallMass_Mass = 1000;
        reco_WSmallMass_CosJets = 1000;
        reco_Cos = 1000;
    }
    
    
          
    tree->Fill();

    nEventsProcessed++;
    if (nEventsProcessed % 10000 == 0) {
        std::cout << nEventsProcessed << " events processed" << std::endl;
    }
}

void MyProcessor::end() {
    tree->Write();
    file->Close();
}
