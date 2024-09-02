#ifndef DELPHES_TOOLS_GENPARTICLE_H_
#define DELPHES_TOOLS_GENPARTICLE_H_


#include "classes/DelphesClasses.h"
#include <vector>


using namespace ROOT::Math;


std::vector<const GenParticle*> getDaughters(const GenParticle* mother, const std::vector<const GenParticle*>& p_arr);

const GenParticle* getHardProcessMotherCopy(const GenParticle*,
                                            const std::vector<const GenParticle*>&);

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L511-L523
const GenParticle* getPreviousCopy(const GenParticle* p,
                                   const std::vector<const GenParticle*>& p_arr);

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L525-L537
const GenParticle* getNextCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_arr);

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L392-L404
const GenParticle* getFirstCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_arr);

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L406-L418
const GenParticle* getLastCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_arr);

bool isFirstCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_arr);

bool isLastCopy(const GenParticle* p, const std::vector<const GenParticle*>& p_arr);

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L265-L295
// TODO
bool isHardProcess(const GenParticle* p);

bool isFromHardProcess(const GenParticle*, const std::vector<const GenParticle*>&);

const GenParticle* findDaughter(const GenParticle* mother,
                                const int pid,
                                const std::vector<const GenParticle*>& p_arr);

const GenParticle* findLastDaughter(const GenParticle* mother,
                                const int pid,
                                const std::vector<const GenParticle*>& p_arr);


// match parton to jet with shortest distance
// starting with the shortest distance available
// apply some outlier rejection if desired
// https://github.com/cms-sw/cmssw/blob/CMSSW_14_1_0_pre4/TopQuarkAnalysis/TopTools/src/JetPartonMatching.cc#L149-L196
std::vector<int32_t> matchPartonToJet(const std::vector<const GenParticle*> parton_vec,
                                        const std::vector<const Jet*> jet_vec,
                                        const float max_distance);





#endif // DELPHES_TOOLS_GENPARTICLE_H_
