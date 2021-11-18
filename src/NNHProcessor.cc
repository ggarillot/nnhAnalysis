#include "NNHProcessor.hh"

#include "EventShape.hh"
#include "ParticleInfo.hh"

#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <marlin/VerbosityLevels.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

#include <TFile.h>
#include <TTree.h>

#include <fastjet/ClusterSequence.hh>

#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/ThreeVector.h>

#include <streamlog/loglevels.h>
#include <streamlog/streamlog.h>
#include <string>
#include <vector>

NNHProcessor aNNHProcessor;

double dijDistance(const fastjet::PseudoJet& i, const fastjet::PseudoJet& j)
{
    CLHEP::Hep3Vector p_i(i.px(), i.py(), i.pz());
    CLHEP::Hep3Vector p_j(j.px(), j.py(), j.pz());

    auto cos_ij = p_i.dot(p_j) / (p_i.mag() * p_j.mag());
    auto minE = std::min(i.e(), j.e());

    return 2.0 * minE * (1 - cos_ij);
}

double yijDistance(const fastjet::PseudoJet& i, const fastjet::PseudoJet& j, const double Q)
{
    return dijDistance(i, j) / (Q * Q);
}

NNHProcessor::NNHProcessor()
    : Processor("NNHProcessor")
{
    registerProcessorParameter("RootFileName", "File name for the root output", rootFileName, std::string("test.root"));

    registerProcessorParameter("MCParticlesCollectionName", "Name of the MC particles collection",
                               mcParticleCollectionName, std::string("MCParticlesSkimmed"));

    registerProcessorParameter("ReconstructedParticlesCollectionName", "Name of the reconstructed particles collection",
                               reconstructedParticleCollectionName, std::string("PandoraPFOs"));

    registerProcessorParameter("IsolatedLeptonsCollectionNames",
                               "Name of the reconstructed isolated leptons collections", isolatedLeptonsCollectionNames,
                               {"IsolatedElectrons", "IsolatedMuons", "IsolatedTaus"});
}

void NNHProcessor::init()
{
    streamlog_out(MESSAGE) << "NNHProcessor::init()" << std::endl;
    outputFile = new TFile(rootFileName.c_str(), "RECREATE");
    outputTree = new TTree("tree", "tree");

    // Collision Energy
    outputTree->Branch("processID", &processID);
    outputTree->Branch("event", &event);
    outputTree->Branch("sqrtS", &sqrtS);

    outputTree->Branch("isValid", &isValid);
    outputTree->Branch("visible_e", &visible_e);
    outputTree->Branch("visible_pt", &visible_pt);
    outputTree->Branch("visible_m", &visible_m);
    outputTree->Branch("visible_recMass", &visible_recMass);
    outputTree->Branch("nParticles", &nParticles);

    outputTree->Branch("nIsoLep", &nIsoLep);

    outputTree->Branch("w1_m", &w1_m);
    outputTree->Branch("w1_pt", &w1_pt);
    outputTree->Branch("w1_e", &w1_e);
    outputTree->Branch("w1_cosBetw", &w1_cosBetw);
    outputTree->Branch("w2_m", &w2_m);
    outputTree->Branch("w2_pt", &w2_pt);
    outputTree->Branch("w2_e", &w2_e);
    outputTree->Branch("w2_cosBetw", &w2_cosBetw);

    outputTree->Branch("higgs_cosBetw", &higgs_cosBetw);

    outputTree->Branch("y_12", &y_12);
    outputTree->Branch("y_23", &y_23);
    outputTree->Branch("y_34", &y_34);
    outputTree->Branch("y_45", &y_45);
    outputTree->Branch("y_56", &y_56);
    outputTree->Branch("y_67", &y_67);

    outputTree->Branch("zz_z1_m", &zz_z1_m);
    outputTree->Branch("zz_z2_m", &zz_z2_m);
    outputTree->Branch("sl_w_m", &sl_w_m);
    outputTree->Branch("sl_rec_m", &sl_rec_m);

    outputTree->Branch("sphericity", &sphericity);
    outputTree->Branch("cosThrust", &cosThrust);
    outputTree->Branch("majorThrust", &majorThrust);
    outputTree->Branch("minorThrust", &minorThrust);

    // ISR
    outputTree->Branch("mc_ISR_e", &mc_ISR_e);
    outputTree->Branch("mc_ISR_pt", &mc_ISR_pt);

    // Neutrinos
    outputTree->Branch("mc_nu_flavor", &mc_nu_flavor);
    outputTree->Branch("mc_nu_e", &mc_nu_e);
    outputTree->Branch("mc_nu_pt", &mc_nu_pt);
    outputTree->Branch("mc_nu_m", &mc_nu_m);
    outputTree->Branch("mc_nu_cosBetw", &mc_nu_cosBetw);

    // Higgs
    outputTree->Branch("mc_higgs_e", &mc_higgs_e);
    outputTree->Branch("mc_higgs_pt", &mc_higgs_pt);
    outputTree->Branch("mc_higgs_m", &mc_higgs_m);
    outputTree->Branch("mc_higgs_recMass", &mc_higgs_recMass);
    outputTree->Branch("mc_higgs_decay", &mc_higgs_decay);
    outputTree->Branch("mc_higgs_subDecay", &mc_higgs_subDecay);

    outputTree->Branch("mc_higgs_decay1_e", &mc_higgs_decay1_e);
    outputTree->Branch("mc_higgs_decay1_pt", &mc_higgs_decay1_pt);
    outputTree->Branch("mc_higgs_decay1_m", &mc_higgs_decay1_m);
    outputTree->Branch("mc_higgs_decay2_e", &mc_higgs_decay2_e);
    outputTree->Branch("mc_higgs_decay2_pt", &mc_higgs_decay2_pt);
    outputTree->Branch("mc_higgs_decay2_m", &mc_higgs_decay2_m);
    outputTree->Branch("mc_higgs_decay_cosBetw", &mc_higgs_decay_cosBetw);

    // // Energy and Pt of the collision
    // outputTree->Branch("mc_TotalEnergy", &mc_TotalEnergy);
    // outputTree->Branch("mc_TotalPT", &mc_TotalPT);

    // // Number of jets
    // outputTree->Branch("reco_njets", &reco_njets);

    // // Reconstructed jets parameters
    // outputTree->Branch("reco_Y12", &reco_Y12);
    // outputTree->Branch("reco_Y13", &reco_Y13);
    // outputTree->Branch("reco_Y14", &reco_Y14);
    // outputTree->Branch("reco_Y23", &reco_Y23);
    // outputTree->Branch("reco_Y24", &reco_Y24);
    // outputTree->Branch("reco_Y34", &reco_Y34);
}

void NNHProcessor::clear()
{
    isValid = false;
    particles.clear();
    isolatedLeptons.clear();
}

fastjet::PseudoJet recoParticleToPseudoJet(EVENT::ReconstructedParticle* recoPart)
{
    auto mom = recoPart->getMomentum();
    auto energy = recoPart->getEnergy();

    fastjet::PseudoJet particle(mom[0], mom[1], mom[2], energy);
    auto               partInfo = new ParticleInfo;
    partInfo->setRecoParticle(recoPart);
    particle.set_user_info(partInfo);
    return particle;
}

double computeRecoilMass(const CLHEP::HepLorentzVector z4Vector, float energy)
{
    auto pTot = CLHEP::Hep3Vector(energy * std::sin(7e-3), 0, 0);
    pTot = pTot - CLHEP::Hep3Vector(z4Vector.px(), z4Vector.py(), z4Vector.pz());
    double rm = (energy - z4Vector.e()) * (energy - z4Vector.e()) - pTot.mag2();

    if (rm < 0)
        throw(std::logic_error("Impossible recoil mass : m2<0"));

    rm = std::sqrt(rm);
    return rm;
}
double computeRecoilMass(const fastjet::PseudoJet& particle, float energy)
{
    const auto vec = CLHEP::HepLorentzVector(particle.px(), particle.py(), particle.pz(), particle.e());
    return computeRecoilMass(vec, energy);
}

void NNHProcessor::processISR(const EVENT::MCParticle* gamma0, const EVENT::MCParticle* gamma1)
{
    mc_ISR_e = -1;
    mc_ISR_pt = -1;

    if (gamma0->getPDG() != 22 || gamma1->getPDG() != 22)
        throw std::logic_error("not gammas");

    auto gamma0_4Vec = CLHEP::HepLorentzVector(gamma0->getMomentum()[0], gamma0->getMomentum()[1],
                                               gamma0->getMomentum()[2], gamma0->getEnergy());
    auto gamma1_4Vec = CLHEP::HepLorentzVector(gamma1->getMomentum()[0], gamma1->getMomentum()[1],
                                               gamma1->getMomentum()[2], gamma1->getEnergy());
    auto ISR_4Vec = gamma0_4Vec + gamma1_4Vec;

    mc_ISR_e = ISR_4Vec.e();
    mc_ISR_pt = ISR_4Vec.perp();
}

void NNHProcessor::processNeutrinos(const EVENT::MCParticle* nu0, const EVENT::MCParticle* nu1)
{
    mc_nu_flavor = -1;
    mc_nu_e = -1;
    mc_nu_pt = -1;
    mc_nu_m = -1;
    mc_nu_cosBetw = -2;

    auto nu0_flavor = nu0->getPDG();
    auto nu1_flavor = nu1->getPDG();

    auto flavor = std::abs(nu0_flavor);

    bool isNeutrinos = (nu0_flavor == -nu1_flavor) && (flavor == 12 || flavor == 14 || flavor == 16);

    if (!isNeutrinos)
        throw std::logic_error("not neutrinos");

    auto nu0_4Vec =
        CLHEP::HepLorentzVector(nu0->getMomentum()[0], nu0->getMomentum()[1], nu0->getMomentum()[2], nu0->getEnergy());
    auto nu1_4Vec =
        CLHEP::HepLorentzVector(nu1->getMomentum()[0], nu1->getMomentum()[1], nu1->getMomentum()[2], nu1->getEnergy());
    auto angleBetw = nu0_4Vec.v().angle(nu1_4Vec.v());

    mc_nu_flavor = flavor;
    mc_nu_e = (nu0_4Vec + nu1_4Vec).e();
    mc_nu_pt = (nu0_4Vec + nu1_4Vec).perp();
    mc_nu_m = (nu0_4Vec + nu1_4Vec).m();
    mc_nu_cosBetw = std::cos(angleBetw);
}

void NNHProcessor::processHiggs(const EVENT::MCParticle* higgs)
{
    mc_higgs_e = -1;
    mc_higgs_pt = -1;
    mc_higgs_m = -1;
    mc_higgs_recMass = -1;
    mc_higgs_decay = -1;
    mc_higgs_subDecay = -1;

    mc_higgs_decay1_e = -1;
    mc_higgs_decay1_pt = -1;
    mc_higgs_decay1_m = -1;
    mc_higgs_decay2_e = -1;
    mc_higgs_decay2_pt = -1;
    mc_higgs_decay2_m = -1;
    mc_higgs_decay_cosBetw = -1;

    if (higgs->getPDG() != 25)
        throw std::logic_error("not a higgs");

    auto higgs_4Vec = CLHEP::HepLorentzVector(higgs->getMomentum()[0], higgs->getMomentum()[1], higgs->getMomentum()[2],
                                              higgs->getEnergy());

    auto vec = higgs->getDaughters();

    if (vec.size() != 2)
        throw std::logic_error("weird higgs decay : not 2 particles");

    auto decay = findDecayMode(vec[0], vec[1]);

    auto decay1_4Vec = CLHEP::HepLorentzVector(vec[0]->getMomentum()[0], vec[0]->getMomentum()[1],
                                               vec[0]->getMomentum()[2], vec[0]->getEnergy());
    auto decay2_4Vec = CLHEP::HepLorentzVector(vec[1]->getMomentum()[0], vec[1]->getMomentum()[1],
                                               vec[1]->getMomentum()[2], vec[1]->getEnergy());

    if (decay2_4Vec.e() > decay1_4Vec.e())
        std::swap(decay1_4Vec, decay2_4Vec);

    auto angleBetw = decay1_4Vec.v().angle(decay2_4Vec.v());

    mc_higgs_e = higgs_4Vec.e();
    mc_higgs_pt = higgs_4Vec.perp();
    mc_higgs_m = higgs_4Vec.m();

    try
    {
        mc_higgs_recMass = computeRecoilMass(higgs_4Vec, sqrtS);
    }
    catch (std::logic_error& e)
    {
        mc_higgs_recMass = 0;
    }

    mc_higgs_decay = decay[0];
    mc_higgs_subDecay = decay[1];

    mc_higgs_decay1_e = decay1_4Vec.e();
    mc_higgs_decay1_pt = decay1_4Vec.perp();
    mc_higgs_decay1_m = decay1_4Vec.m();
    mc_higgs_decay2_e = decay2_4Vec.e();
    mc_higgs_decay2_pt = decay2_4Vec.perp();
    mc_higgs_decay2_m = decay2_4Vec.m();
    mc_higgs_decay_cosBetw = std::cos(angleBetw);
}

std::array<fastjet::PseudoJet, 2> findParticleByMass(const std::vector<fastjet::PseudoJet> jets,
                                                     const double                          targetMass,
                                                     std::vector<fastjet::PseudoJet>&      remainingJets)
{
    auto goodPair = std::array<unsigned int, 2>{};
    auto chi2 = std::numeric_limits<double>::max();

    for (auto i = 0U; i < jets.size(); ++i)
    {
        for (auto j = i + 1; j < jets.size(); ++j)
        {
            auto m = (jets[i] + jets[j]).m();
            auto val = std::abs(m - targetMass);

            if (val > chi2)
                continue;

            chi2 = val;
            goodPair = {i, j};
        }
    }

    std::array<fastjet::PseudoJet, 2> toReturn = {jets[goodPair[0]], jets[goodPair[1]]};

    for (auto i = 0U; i < jets.size(); ++i)
    {
        if (i == goodPair[0] || i == goodPair[1])
            continue;
        remainingJets.push_back(jets[i]);
    }

    return toReturn;
}

std::array<int, 2> NNHProcessor::findDecayMode(const EVENT::MCParticle* part1, const EVENT::MCParticle* part2) const
{
    std::array<int, 2> toReturn{-1, -1};

    auto decay1 = std::abs(part1->getPDG());
    auto decay2 = std::abs(part2->getPDG());

    if (decay1 != 22 && decay1 != 23)
    {
        if (decay1 != decay2)
            throw std::logic_error("weird higgs decay : " + std::to_string(decay1) + ", " + std::to_string(decay2));
    }

    if (decay1 != decay2)
        decay1 = 25;

    decay2 = 0;
    if (decay1 == 24 || decay1 == 23)
    {
        auto vec0 = part1->getDaughters();
        auto vec1 = part2->getDaughters();

        if (vec0.size() != 2 || vec1.size() != 2)
            throw std::logic_error("weird higgs subdecay for WW or ZZ : " + std::to_string(vec0.size() + vec1.size()) +
                                   "particles");

        std::array<int, 4> subDecay = {{std::abs(vec0[0]->getPDG()), std::abs(vec0[1]->getPDG()),
                                        std::abs(vec1[0]->getPDG()), std::abs(vec1[1]->getPDG())}};

        std::sort(subDecay.begin(), subDecay.end());

        if (subDecay[1] < 10) // qq--
        {
            if (subDecay[3] < 10) // qqqq
                decay2 = 1;
            else if (subDecay[2] % 2 != 0 && subDecay[3] % 2 != 0) // qqll
            {
                if (subDecay[2] == 11)
                    decay2 = 21;
                else if (subDecay[2] == 13)
                    decay2 = 22;
                else if (subDecay[2] == 15)
                    decay2 = 23;
                else
                    throw std::logic_error("weird qqll decay");
            }
            else if (subDecay[2] % 2 == 0 && subDecay[3] % 2 == 0) // qqvv
            {
                decay2 = 4;
            }
            else // qqlv
            {
                if (subDecay[2] == 11)
                    decay2 = 31;
                else if (subDecay[2] == 13)
                    decay2 = 32;
                else if (subDecay[2] == 15)
                    decay2 = 33;
                else
                    throw std::logic_error("weird qqlv decay");
            }
        }
        else
        {
            int nbNu = 0;
            for (const auto& i : subDecay)
            {
                if (i % 2 == 0)
                    nbNu++;
            }

            if (nbNu == 0) // llll
            {
                decay2 = 500;
                if (subDecay[0] == 11)
                    decay2 += 10;
                else if (subDecay[0] == 13)
                    decay2 += 20;
                else if (subDecay[0] == 15)
                    decay2 += 30;
                else
                    throw std::logic_error("weird llll decay");

                if (subDecay[2] == 11)
                    decay2 += 1;
                else if (subDecay[2] == 13)
                    decay2 += 2;
                else if (subDecay[2] == 15)
                    decay2 += 3;
                else
                    throw std::logic_error("weird llll decay");
            }
            else if (nbNu == 2) // llvv
            {
                decay2 = 600;
                std::vector<int> temp = {};
                for (const auto& i : subDecay)
                {
                    if (i % 2 != 0)
                    {
                        if (i == 11)
                            temp.push_back(1);
                        else if (i == 13)
                            temp.push_back(2);
                        else if (i == 15)
                            temp.push_back(3);
                        else
                            throw std::logic_error("weird llvv decay");
                    }
                }

                if (temp.size() != 2)
                    throw std::logic_error("weird llvv decay");

                std::sort(temp.begin(), temp.end());

                decay2 = decay2 + 10 * temp[0] + temp[1];
            }
            else // vvvv
            {
                decay2 = 7;
            }
        }
    }

    toReturn = {decay1, decay2};
    return toReturn;
}

Eigen::Matrix3d NNHProcessor::computeSphericityTensor(const std::vector<fastjet::PseudoJet>& particleVec) const
{
    Eigen::Matrix3d tensor;

    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            double num = 0;
            double denom = 0;

            for (const auto& particle : particleVec)
            {
                num += particle.four_mom()[i] * particle.four_mom()[j];
                denom += particle.modp2();
            }
            tensor(i, j) = num / denom;
        }
    }

    return tensor;
}

double NNHProcessor::computeSphericity(const std::vector<fastjet::PseudoJet>& particleVec) const
{
    auto tensor = computeSphericityTensor(particleVec);

    auto eigenVal = tensor.eigenvalues();

    std::array<double, 3> val = {{std::norm(eigenVal(0)), std::norm(eigenVal(1)), std::norm(eigenVal(2))}};
    std::sort(val.begin(), val.end());

    streamlog_out(DEBUG) << "Sphericity eigenvalues : (" << val[0] << " " << val[1] << " " << val[2] << ")"
                         << std::endl;

    return 1.5 * (val[0] + val[1]);
}

void NNHProcessor::processEvent(LCEvent* evt)
{
    clear();

    processID = evt->getParameters().getIntVal(std::string("ProcessID"));
    event = evt->getParameters().getIntVal(std::string("Event Number"));
    sqrtS = evt->getParameters().getFloatVal(std::string("Energy"));

    for (const auto& colName : isolatedLeptonsCollectionNames)
    {
        auto col = evt->getCollection(colName);
        auto n = col->getNumberOfElements();

        for (auto i = 0; i < n; ++i)
        {
            auto particle = dynamic_cast<EVENT::ReconstructedParticle*>(col->getElementAt(i));
            isolatedLeptons.insert(particle);
        }
    }

    nIsoLep = isolatedLeptons.size();

    recoCol = evt->getCollection(reconstructedParticleCollectionName);
    nParticles = recoCol->getNumberOfElements();
    particles.reserve(nParticles);

    for (int index = 0; index < nParticles; ++index)
    {
        auto recoPart = dynamic_cast<EVENT::ReconstructedParticle*>(recoCol->getElementAt(index));

        auto particle = recoParticleToPseudoJet(recoPart);
        particles.push_back(particle);
    }

    // MC stuff
    mcCol = evt->getCollection(mcParticleCollectionName);

    auto gamma0 = dynamic_cast<EVENT::MCParticle*>(mcCol->getElementAt(6));
    auto gamma1 = dynamic_cast<EVENT::MCParticle*>(mcCol->getElementAt(7));
    auto nu0 = dynamic_cast<EVENT::MCParticle*>(mcCol->getElementAt(8));
    auto nu1 = dynamic_cast<EVENT::MCParticle*>(mcCol->getElementAt(9));
    auto higgs = dynamic_cast<EVENT::MCParticle*>(mcCol->getElementAt(10));

    try
    {
        processISR(gamma0, gamma1);
    }
    catch (std::logic_error& e)
    {
        streamlog_out(DEBUG) << "Run : " << evt->getRunNumber() << ", Event : " << evt->getEventNumber() << " : "
                             << e.what() << std::endl;
    }
    try
    {
        processNeutrinos(nu0, nu1);
    }
    catch (std::logic_error& e)
    {
        streamlog_out(DEBUG) << "Run : " << evt->getRunNumber() << ", Event : " << evt->getEventNumber() << " : "
                             << e.what() << std::endl;
    }
    try
    {
        processHiggs(higgs);
    }
    catch (std::logic_error& e)
    {
        streamlog_out(DEBUG) << "Run : " << evt->getRunNumber() << ", Event : " << evt->getEventNumber() << " : "
                             << e.what() << std::endl;
    }

    // end of MC stuff

    EventShape shape;
    shape.setPartList(particles);

    cosThrust = std::abs(shape.thrustAxis().cosTheta());
    majorThrust = shape.majorThrust();
    minorThrust = shape.minorThrust();

    if (minorThrust != minorThrust) // handle NaN case
        minorThrust = 0;

    sphericity = computeSphericity(particles);

    // Jets study

    std::array<std::vector<fastjet::PseudoJet>, 7> jets_perN{};
    std::array<float, 7>                           yCutArray{};

    auto j = fastjet::JetDefinition(fastjet::ee_kt_algorithm);
    auto cs = fastjet::ClusterSequence(particles, j);

    auto maxJets = std::min(6, nParticles);

    if (maxJets == 0)
    {
        streamlog_out(WARNING) << "Run : " << evt->getRunNumber() << ", Event : " << evt->getEventNumber() << " : "
                               << "No particles in event !" << std::endl;
    }

    for (int i = 1; i <= maxJets; ++i)
    {
        auto jets = sorted_by_E(cs.exclusive_jets(i));
        jets_perN[i] = jets;

        yCutArray[i] = cs.exclusive_ymerge(i);
    }

    //---------------------

    // std::vector<fastjet::PseudoJet> jets_final = cs.exclusive_jets_ycut(y_cut);
    // reco_njets = jets_final.size();

    if (!jets_perN[1].empty())
    {
        visible_e = jets_perN[1][0].e();
        visible_pt = jets_perN[1][0].pt();
        visible_m = jets_perN[1][0].m();
        try
        {
            visible_recMass = computeRecoilMass(jets_perN[1][0], sqrtS);
        }
        catch (std::logic_error& e)
        {
            visible_recMass = 0;
        }
    }

    // 3 jets study
    if (!jets_perN[3].empty())
    {
        const auto& jets = jets_perN[3];

        std::vector<fastjet::PseudoJet> osef{};

        auto W_jetPair = findParticleByMass(jets, W_MASS_REF, osef);
        auto W = join(W_jetPair[0], W_jetPair[1]);

        sl_w_m = W.m();

        try
        {
            sl_rec_m = computeRecoilMass(W, sqrtS);
        }
        catch (std::logic_error& e)
        {
            streamlog_out(DEBUG) << "Run : " << evt->getRunNumber() << ", Event : " << evt->getEventNumber() << " : "
                                 << "Impossible recoil mass !" << std::endl;
            sl_rec_m = 0;
        }
    }

    // 4 jets study
    if (!jets_perN[4].empty())
    {
        isValid = true;

        const auto& jets = jets_perN[4];

        std::vector<fastjet::PseudoJet> smallW_jetPair{};

        auto bigW_jetPair = findParticleByMass(jets, W_MASS_REF, smallW_jetPair);

        auto bigW = join(bigW_jetPair[0], bigW_jetPair[1]);
        auto bigW_mom = CLHEP::Hep3Vector(bigW.px(), bigW.py(), bigW.pz());
        auto bigW_jet1Mom = CLHEP::Hep3Vector(bigW_jetPair[0].px(), bigW_jetPair[0].py(), bigW_jetPair[0].pz());
        auto bigW_jet2Mom = CLHEP::Hep3Vector(bigW_jetPair[1].px(), bigW_jetPair[1].py(), bigW_jetPair[1].pz());

        auto smallW = join(smallW_jetPair[0], smallW_jetPair[1]);
        auto smallW_mom = CLHEP::Hep3Vector(smallW.px(), smallW.py(), smallW.pz());
        auto smallW_jet1Mom = CLHEP::Hep3Vector(smallW_jetPair[0].px(), smallW_jetPair[0].py(), smallW_jetPair[0].pz());
        auto smallW_jet2Mom = CLHEP::Hep3Vector(smallW_jetPair[1].px(), smallW_jetPair[1].py(), smallW_jetPair[1].pz());

        w1_m = bigW.m();
        w1_pt = bigW.pt();
        w1_e = bigW.e();
        w1_cosBetw = std::cos(bigW_jet1Mom.angle(bigW_jet2Mom));

        w2_m = smallW.m();
        w2_pt = smallW.pt();
        w2_e = smallW.e();
        w2_cosBetw = std::cos(smallW_jet1Mom.angle(smallW_jet2Mom));

        higgs_cosBetw = std::cos(bigW_mom.angle(smallW_mom));

        // background study
        std::vector<fastjet::PseudoJet> smallZ_jetPair{};

        auto bigZ_jetPair = findParticleByMass(jets, Z_MASS_REF, smallZ_jetPair);

        auto bigZ = bigZ_jetPair[0] + bigZ_jetPair[1];
        auto smallZ = smallZ_jetPair[0] + smallZ_jetPair[1];

        zz_z1_m = bigZ.m();
        zz_z2_m = smallZ.m();
    }

    // if (nRecoParticles > 3)
    // { // we need at least 4 reconstructed particles to study Y_ij distributions
    //     reco_Y12 = jet_parameter(jets_final_bis[0], jets_final_bis[1], 250);
    //     reco_Y13 = jet_parameter(jets_final_bis[0], jets_final_bis[2], 250);
    //     reco_Y14 = jet_parameter(jets_final_bis[0], jets_final_bis[3], 250);
    //     reco_Y23 = jet_parameter(jets_final_bis[1], jets_final_bis[2], 250);
    //     reco_Y24 = jet_parameter(jets_final_bis[1], jets_final_bis[3], 250);
    //     reco_Y34 = jet_parameter(jets_final_bis[2], jets_final_bis[3], 250);

    constexpr float minYCut = std::numeric_limits<float>::min();
    for (auto& yCut : yCutArray)
        yCut = std::max(yCut, minYCut);

    y_12 = -log10(yCutArray[1]);
    y_23 = -log10(yCutArray[2]);
    y_34 = -log10(yCutArray[3]);
    y_45 = -log10(yCutArray[4]);
    y_56 = -log10(yCutArray[5]);
    y_67 = -log10(yCutArray[6]);

    outputTree->Fill();

    nEventsProcessed++;
    if (nEventsProcessed % 10000 == 0)
        streamlog_out(MESSAGE) << nEventsProcessed << " events processed" << std::endl;
}

void NNHProcessor::end()
{
    outputTree->Write();
    outputFile->Close();
}