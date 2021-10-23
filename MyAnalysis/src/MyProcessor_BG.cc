#include "MyProcessor_BG.h"

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

MyProcessor_BG aMyProcessor_BG;

double jet_parameter_1(const fastjet::PseudoJet& jet_1, const fastjet::PseudoJet& jet_2, double total_E) {   
    CLHEP::Hep3Vector p1(jet_1.px(), jet_1.py(), jet_1.pz());
    CLHEP::Hep3Vector p2(jet_2.px(), jet_2.py(), jet_2.pz());
    double log_Y = (-1)*std::log(jet_1.E()*jet_2.E() * (1 - (p1.dot(p2))/(p1.mag()*p2.mag()))* 2./(total_E*total_E));  
    return log_Y;
}

double dist_mass_1(const fastjet::PseudoJet& jet_1, const fastjet::PseudoJet& jet_2) {
    double d = std::min(std::abs(jet_1.m() - 80.379) + std::abs(jet_2.m() - 37.654), std::abs(jet_2.m() - 80.379) + std::abs(jet_1.m() - 37.654));
    return d;
}

bool test_order_1(const fastjet::PseudoJet& jet_1, const fastjet::PseudoJet& jet_2, const double d) {
    if (d == std::abs(jet_1.m() - 80.379) + std::abs(jet_2.m() - 37.654)) {
        bool order = true;
        return order;
    } else {
        bool order = false;
        return order;
    }
}

double cos_1(const fastjet::PseudoJet& jet_1, const fastjet::PseudoJet& jet_2) {
    double cos = (jet_1.px()*jet_2.px() + jet_1.py()*jet_2.py() + jet_1.pz()*jet_2.pz())/(jet_1.modp()*jet_2.modp());
    return cos;
}

MyProcessor_BG::MyProcessor_BG()
    : Processor("MyProcessor_BG")
{
    registerProcessorParameter("RootFileName", "File name for the root output", rootFileName, std::string("test.root"));

    registerProcessorParameter("MCParticlesCollectionName", "Name of the MC particles collection", mcParticleCollectionName, std::string("MCParticlesSkimmed"));
    
    registerProcessorParameter("ReconstructedParticlesCollectionName", "Name of the reconstructed particles collection",reconstructedParticleCollectionName, std::string("PandoraPFOs"));
}

void MyProcessor_BG::init() {
    std::cout << "MyProcessor_BG::init()" << std::endl;
    file = new TFile(rootFileName.c_str(), "RECREATE");
    tree = new TTree("tree", "tree");
    
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

void MyProcessor_BG::clear() {
    reconstructedParticles.clear();
}

void MyProcessor_BG::processEvent(LCEvent* evt) {
    clear();

    //Jets study
    LCCollection* reconstructedParticleCollection = evt->getCollection(reconstructedParticleCollectionName);
    const int     nRecoParticles = reconstructedParticleCollection->getNumberOfElements();
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
        reco_Y12 = jet_parameter_1(jets_final_bis[0], jets_final_bis[1], 250);
        reco_Y13 = jet_parameter_1(jets_final_bis[0], jets_final_bis[2], 250);
        reco_Y14 = jet_parameter_1(jets_final_bis[0], jets_final_bis[3], 250);
        reco_Y23 = jet_parameter_1(jets_final_bis[1], jets_final_bis[2], 250);
        reco_Y24 = jet_parameter_1(jets_final_bis[1], jets_final_bis[3], 250);
        reco_Y34 = jet_parameter_1(jets_final_bis[2], jets_final_bis[3], 250);
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
        std::vector<double> cos_jets_1 = {cos_1(jets_final_bis[0], jets_final_bis[1]), cos_1(jets_final_bis[2], jets_final_bis[3])};
        std::vector<double> cos_jets_2 = {cos_1(jets_final_bis[0], jets_final_bis[2]), cos_1(jets_final_bis[1], jets_final_bis[3])};
        std::vector<double> cos_jets_3 = {cos_1(jets_final_bis[0], jets_final_bis[3]), cos_1(jets_final_bis[1], jets_final_bis[2])};
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
            if (dist_mass_1(jets[i][0], jets[i][2]) < d) {
                d = dist_mass_1(jets[i][0], jets[i][2]);
                order = test_order_1(jets[i][0], jets[i][2], d);
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
        reco_Cos = cos_1(selected_jets[0], selected_jets[1]);
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

void MyProcessor_BG::end() {
    tree->Write();
    file->Close();
}
