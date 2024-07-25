#ifndef DELPHES_TOOLS_SELECTION_H_
#define DELPHES_TOOLS_SELECTION_H_

#include "TClonesArray.h"

#include <vector>
#include <functional>


// Returns selected veto electrons
template<typename T>
std::vector<const T*> runObjectSelection(const TClonesArray* arr,
                                         std::function<bool(const T*)> selection
) {
  std::vector<const T*> output{};
  output.reserve(arr->GetEntries());

  for (int32_t idx = 0; idx < arr->GetEntries(); ++idx) {
    if (auto obj = dynamic_cast<const T*>(arr->At(idx))) {
      if (selection(obj)) {
        output.push_back(obj);
      }

    } else {
      // TODO
      throw std::runtime_error("failed to get");

    }
  }

  return output;
}


#endif // DELPHES_TOOLS_SELECTION_H_
