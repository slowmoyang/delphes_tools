#include "delphes_tools/utils.h"

// ROOT
#include "TSystem.h"

// std
#include <cstdlib>
#include <filesystem>

namespace fs = std::filesystem;


void setupDelphes() {
  if (not std::getenv("DELPHES_PREFIX")) {
    throw std::runtime_error("DELPHES_PREFIX not defined");
  }

  const fs::path delphes_prefix{std::getenv("DELPHES_PREFIX")};
  const fs::path delphes_shared_object = delphes_prefix / "libDelphes.so";
  if (not fs::exists(delphes_shared_object)) {
    throw std::runtime_error("delphes shared object not found: " + delphes_shared_object.string());
  }

  gInterpreter->AddIncludePath(delphes_prefix.c_str());
  gSystem->Load(delphes_shared_object.c_str());
  gInterpreter->Declare("#include \"classes/DelphesClasses.h\"");
}
