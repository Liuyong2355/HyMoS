#ifndef HYMOS_COMMON_SOLVERTYPE_H
#define HYMOS_COMMON_SOLVERTYPE_H

#include <cctype>
#include <string>

#include "ConfigFile.h"

namespace hymos {
namespace common {

enum RuntimeSolverType {
  kSolverSingle = 0,
  kSolverNMG = 1,
  kSolverFMG = 2
};

inline std::string NormalizeSolverToken(const std::string& value) {
  std::string normalized;
  normalized.reserve(value.size());
  for (size_t i = 0; i < value.size(); ++i) {
    const unsigned char ch = static_cast<unsigned char>(value[i]);
    if (std::isalnum(ch)) {
      normalized.push_back(static_cast<char>(std::tolower(ch)));
    }
  }
  return normalized;
}

inline RuntimeSolverType ParseRuntimeSolverType(const ConfigFile& cf,
                                                bool allow_fmg,
                                                RuntimeSolverType default_type = kSolverSingle) {
  std::string raw;
  try {
    raw = (std::string)cf.Value("Solver", "SolverType");
  } catch (const std::string&) {
    return default_type;
  }

  const std::string normalized = NormalizeSolverToken(raw);
  if (normalized == "single") {
    return kSolverSingle;
  }
  if (normalized == "nmg" || normalized == "multigrid") {
    return kSolverNMG;
  }
  if (normalized == "fmg") {
    throw std::string("Solver/SolverType=fmg is not in the current supported scope");
  }
  throw std::string("Unsupported Solver/SolverType = " + raw);
}

inline const char* RuntimeSolverTypeName(RuntimeSolverType type) {
  switch (type) {
    case kSolverNMG:
      return "nmg";
    case kSolverFMG:
      return "fmg";
    case kSolverSingle:
    default:
      return "single";
  }
}

}  // namespace common
}  // namespace hymos

#endif
