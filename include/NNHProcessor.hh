#pragma once

#include <array>
#include <vector>

#include <Eigen/Dense>

#include "EVENT/MCParticle.h"
#include "marlin/Processor.h"

#include <fastjet/PseudoJet.hh>

#include <TFile.h>
#include <TTree.h>

class NNHProcessor : public marlin::Processor
{
  public:
    static constexpr float W_MASS_REF = 80.379;
    static constexpr float Z_MASS_REF = 91.1876;

  public:
    virtual Processor* newProcessor() { return new NNHProcessor; }

    NNHProcessor();
    virtual void init();
    virtual void processEvent(LCEvent* evt);
    virtual void end();

    NNHProcessor(const NNHProcessor& toCopy) = delete;
    void operator=(const NNHProcessor& toCopy) = delete;

  protected:
    void clear();

    void processISR(const EVENT::MCParticle* gamma0, const EVENT::MCParticle* gamma1);
    void processNeutrinos(const EVENT::MCParticle* nu0, const EVENT::MCParticle* nu1);
    void processHiggs(const EVENT::MCParticle* higgs);

    std::array<int, 2> findDecayMode(const EVENT::MCParticle* part1, const EVENT::MCParticle* part2) const;

    Eigen::Matrix3d computeSphericityTensor(const std::vector<fastjet::PseudoJet>& particleVec) const;
    double          computeSphericity(const std::vector<fastjet::PseudoJet>& particleVec) const;

    // steering file parameters
    std::string rootFileName{};
    std::string mcParticleCollectionName{};
    std::string reconstructedParticleCollectionName{};

    TFile* outputFile = nullptr;
    TTree* outputTree = nullptr;

    LCCollection* mcCol = nullptr;
    LCCollection* linkCol = nullptr;
    LCCollection* recoCol = nullptr;

    std::vector<fastjet::PseudoJet> particles{};

    // event variables
    float sqrtS = -1;

    float visible_e = 0;
    float visible_pt = 0;
    float visible_m = 0;

    bool isValid = false;

    // Reco variables
    float w1_m = -1;
    float w1_pt = -1;
    float w1_e = -1;
    float w1_cosBetw = -2;

    float w2_m = -1;
    float w2_pt = -1;
    float w2_e = -1;
    float w2_cosBetw = -2;

    float higgs_cosBetw = -2;

    float y_12 = 0;
    float y_23 = 0;
    float y_34 = 0;
    float y_45 = 0;
    float y_56 = 0;
    float y_67 = 0;

    float zz_z1_m = -1;
    float zz_z2_m = -1;
    float sl_w_m = -1;
    float sl_rec_m = -1;

    float sphericity = -1;
    float cosThrust = -2;
    float majorThrust = -1;
    float minorThrust = -1;

    // MC variables
    float mc_ISR_e = 0;
    float mc_ISR_pt = 0;

    int   mc_nu_flavor = 0;
    float mc_nu_e = 0;
    float mc_nu_pt = 0;
    float mc_nu_m = 0;
    float mc_nu_cosBetw = -2;

    float mc_higgs_e = 0;
    float mc_higgs_pt = 0;
    float mc_higgs_m = 0;

    int mc_higgs_decay = 0;
    int mc_higgs_subDecay = 0;

    float mc_higgs_decay1_e = 0;
    float mc_higgs_decay1_pt = 0;
    float mc_higgs_decay1_m = 0;
    float mc_higgs_decay2_e = 0;
    float mc_higgs_decay2_pt = 0;
    float mc_higgs_decay2_m = 0;
    float mc_higgs_decay_cosBetw = -2;

    // debug variables
    int nEventsProcessed = 0;
};
