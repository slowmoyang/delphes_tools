#include "delphes_tools/weight.h"

#include "classes/DelphesClasses.h"

#include <TClonesArray.h>

#include <stdexcept>


float getWeight(const TClonesArray* weight_arr, const size_t idx) {
  return dynamic_cast<const Weight*>(weight_arr->At(idx))->Weight;
}

const std::vector<float> getWeight(const TClonesArray* weight_arr,
                                   const size_t start,
                                   const size_t stop) {
  if (not weight_arr) {
    throw std::runtime_error("weight_arr is nullptr");
  }

  std::vector<float> weight_vec{};
  weight_vec.reserve(stop - start);
  for (size_t idx = start; idx < stop; ++idx) {
    weight_vec.push_back(getWeight(weight_arr, idx));
  }

  return weight_vec;
}

float getCentralWeight(const TClonesArray* weight_arr) {
  return getWeight(weight_arr, 0);
}

const std::vector<float> getScaleWeight(const TClonesArray* weight_arr) {
  if (not weight_arr) {
    throw std::runtime_error("weight_arr is nullptr");
  }
  return getWeight(weight_arr, 1, 9);
}

const std::vector<float> getPDFWeight(const TClonesArray* weight_arr) {
  if (not weight_arr) {
    throw std::runtime_error("weight_arr is nullptr");
  }
  return getWeight(weight_arr, 9, weight_arr->GetEntries());
}

const std::vector<float> getPSWeight(const TClonesArray* weight_arr) {
    if (not weight_arr) {
    throw std::runtime_error("weight_arr is nullptr");
  }
  return getWeight(weight_arr, 0, weight_arr->GetEntries());
}
