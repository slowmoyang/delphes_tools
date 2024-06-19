#include "delphes_tools/genparticle.h"

#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "Math/GenVector/PtEtaPhiM4D.h"
#include "Math/Vector4Dfwd.h"
#include "Math/GenVector/VectorUtil.h"


using namespace ROOT::Math;

const GenParticle* getDaughter(const TClonesArray* gen_particle_array, size_t idx) {
  return dynamic_cast<const GenParticle*>(gen_particle_array->At(idx));
}

std::vector<const GenParticle*> getDaughterVec(const GenParticle* mother, const TClonesArray* p_arr) {
  std::vector<const GenParticle*> output{};
  output.reserve(mother->D2 - mother->D1 + 1);
  for (size_t idx = mother->D1; idx <= mother->D2; ++idx) {
    const GenParticle* daughter = getDaughter(p_arr, idx);
    output.push_back(daughter);
  }
  return output;
}


const GenParticle* getHardProcessMotherCopy(const GenParticle* p,
                                           const TClonesArray* p_arr
) {
  //is particle itself is hard process particle
  if (isHardProcess(p)) {
    return p;

  }

  //check if any other copies are hard process particles
  const GenParticle* pcopy = p;
  std::unordered_set<const GenParticle*> dup_check;
  while (getPreviousCopy(pcopy, p_arr)) {
    dup_check.insert(pcopy);
    pcopy = getPreviousCopy(pcopy, p_arr);
    if (isHardProcess(pcopy)) {
      return pcopy;
    }
    if (dup_check.count(pcopy)) {
      break;
    }
  }
  return nullptr;
}





// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L511-L523
const GenParticle* getPreviousCopy(const GenParticle* p,
                                   const TClonesArray* p_arr) {
  if (p->M1 == -1) {
    throw std::runtime_error("got an unexpected M1=-1");
  }

  if (p->M2 != -1) {
    throw std::runtime_error("got an unexpected M2=" + std::to_string(p->M2));
  }

  auto mother = dynamic_cast<const GenParticle*>(p_arr->At(p->M1));

  return mother->PID == p->PID ? mother : nullptr;
}

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L525-L537
const GenParticle* getNextCopy(const GenParticle* p, const TClonesArray* p_arr) {
  for (size_t idau = p->D1; idau <= p->D2; ++idau) {
    const GenParticle* dau = getDaughter(p_arr, idau);
    if (dau->PID == p->PID) {
      return dau;
    }
  }
  return nullptr;
}

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L392-L404
/////////////////////////////////////////////////////////////////////////////
const GenParticle* getFirstCopy(const GenParticle* p, const TClonesArray* p_arr) {
  const GenParticle* p_copy = p;

  std::unordered_set<const GenParticle*> dup_check;
  while (getPreviousCopy(p_copy, p_arr)) {
    dup_check.insert(p_copy);
    p_copy = getPreviousCopy(p_copy, p_arr);
    if (dup_check.count(p_copy)) {
      return nullptr;
    }
  }
  return p_copy;
}

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L406-L418
const GenParticle* getLastCopy(const GenParticle* p, const TClonesArray* p_arr) {
  const GenParticle* p_copy = p;
  std::unordered_set<const GenParticle*> dup_check;
  while (getNextCopy(p_copy, p_arr)) {
    dup_check.insert(p_copy);
    p_copy = getNextCopy(p_copy, p_arr);
    if (dup_check.count(p_copy)) {
      return nullptr;
    }
  }
  return p_copy;
}

/////////////////////////////////////////////////////////////////////////////
bool isFirstCopy(const GenParticle* p, const TClonesArray* p_arr) {
  return p == getFirstCopy(p, p_arr);
}

bool isLastCopy(const GenParticle* p, const TClonesArray* p_arr) {
  return p == getLastCopy(p, p_arr);
}

// https://github.com/cms-sw/cmssw/blob/CMSSW_13_2_6_patch2/PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h#L265-L295
// TODO
bool isHardProcess(const GenParticle* p) {
  //hard process codes for pythia8 are 21-29 inclusive (currently 21,22,23,24 are used)
  return ((p->Status > 20) and (p->Status < 30));
}


bool isFromHardProcess(const GenParticle* p, const TClonesArray* p_arr) {
  return getHardProcessMotherCopy(p, p_arr) != nullptr;
}


const GenParticle* findDaughter(const GenParticle* mother,
                                const int pid,
                                const TClonesArray* p_arr) {
  mother = getLastCopy(mother, p_arr);
  for (const GenParticle* daughter : getDaughterVec(mother, p_arr)) {
    if (daughter->PID == pid) {
      return daughter;
    }
  }

  return nullptr;
}

const GenParticle* findLastDaughter(const GenParticle* mother,
                                const int pid,
                                const TClonesArray* p_arr
) {
  const GenParticle* daughter = findDaughter(mother, pid, p_arr);
  if (daughter) {
    daughter = getLastCopy(daughter, p_arr);
  }
  return daughter;
}


// match parton to jet with shortest distance
// starting with the shortest distance available
// apply some outlier rejection if desired
// https://github.com/cms-sw/cmssw/blob/CMSSW_14_1_0_pre4/TopQuarkAnalysis/TopTools/src/JetPartonMatching.cc#L149-L196
std::vector<int32_t> matchPartonToJet(const std::vector<const GenParticle*> parton_vec,
                                        const std::vector<const Jet*> jet_vec,
                                        const float max_distance) {
  // FIXME update comments
  // prepare vector of pairs with distances between
  // all partons to all jets in the input vectors
  std::vector<std::tuple<double, size_t, size_t> > match_vec{};
  match_vec.reserve(parton_vec.size() * jet_vec.size());

  for (size_t parton_idx = 0; parton_idx < parton_vec.size(); ++parton_idx) {
    for (size_t jet_idx = 0; jet_idx < jet_vec.size(); ++jet_idx) {

      const double distance = VectorUtil::DeltaR(
          jet_vec.at(jet_idx)->P4(),
          parton_vec.at(parton_idx)->P4());
      match_vec.emplace_back(distance, parton_idx, jet_idx);
    }
  }
  std::sort(match_vec.begin(), match_vec.end());

  std::set<size_t> seen;

  std::vector<int32_t> parton_jet_idx(parton_vec.size(), -1);

  while (not match_vec.empty()) {
    auto [distance, parton_idx, jet_idx] = match_vec.at(0);
    // use primitive outlier rejection if desired
    if (distance > max_distance) {
      jet_idx = -1;
    }

    if (seen.contains(parton_idx)) {
      throw std::runtime_error("got a parton already seen");
    }

    parton_jet_idx.at(parton_idx) = jet_idx;
    seen.insert(parton_idx);

    // remove all values for the matched parton
    // and the matched jet
    size_t other_parton_idx;
    size_t other_jet_idx;
    // XXX
    // the original code use unsigned int
    for (int other_idx = 0; other_idx < match_vec.size(); ++other_idx) {
      if (other_idx < 0) {
        throw std::runtime_error("got a netgative other_idx: " + std::to_string(other_idx));
      }

      std::tie(std::ignore, other_parton_idx, other_jet_idx) = match_vec.at(other_idx);
      if ((other_parton_idx == parton_idx) or (other_jet_idx == jet_idx)) {
        match_vec.erase(match_vec.begin() + other_idx);
        --other_idx;
      }
    }
  }

  return parton_jet_idx;
}


/////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////



