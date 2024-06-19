#ifndef DELPHES_TOOLS_WEIGHT_H_
#define DELPHES_TOOLS_WEIGHT_H_

#include <cstddef>
#include <vector>

class TClonesArray;

float getWeight(const TClonesArray* weight_arr, const size_t idx);

const std::vector<float> getWeight(const TClonesArray* weight_arr,
                                   const size_t start,
                                   const size_t stop);

float getCentralWeight(const TClonesArray* weight_arr);

const std::vector<float> getScaleWeight(const TClonesArray* weight_arr);

const std::vector<float> getPDFWeight(const TClonesArray* weight_arr);

const std::vector<float> getPSWeight(const TClonesArray* weight_arr);

#endif // DELPHES_TOOLS_WEIGHT_H_
