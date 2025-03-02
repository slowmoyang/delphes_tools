#include "delphes_tools/utils.h"

// ROOT
#include "TSystem.h"

// std
#include <cstdlib>
#include <filesystem>

namespace fs = std::filesystem;


void setupDelphes(const std::string env_name) {
  if (not std::getenv(env_name.c_str())) {
    throw std::runtime_error(env_name + " not defined");
  }

  const fs::path prefix{std::getenv(env_name.c_str())};
  if (not fs::exists(prefix)) {
    throw std::runtime_error("Prefix does not exist: " + prefix.string());
  }

  fs::path so_file, include_dir;
  if (env_name == "DELPHES_PREFIX") {
    so_file = prefix / "libDelphes.so";
    include_dir = prefix;

  } else if (env_name == "CONDA_PREFIX") {
    so_file = prefix / "lib" / "libDelphes.so";
    include_dir = prefix / "include";

  } else {
    throw std::runtime_error("Unknown environment variable: " + env_name);

  }

  if (not fs::exists(so_file)) {
    throw std::runtime_error("shared object not found: " + so_file.string());
  }

  if (not fs::exists(include_dir)) {
    throw std::runtime_error("include directory not found: " + include_dir.string());
  }

  gInterpreter->AddIncludePath(include_dir.c_str());
  gSystem->Load(so_file.c_str());
  gInterpreter->Declare("#include \"classes/DelphesClasses.h\"");
}
