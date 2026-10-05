// Usage:
// - This file is the single place where built-in cases are registered.
// - To add a new backend case, add one `CaseEntry` here and keep the case code
//   itself inside `cases/`.

#include "CaseRegistry.h"

#include <cctype>

#include "BuiltinLaunchers.h"

namespace hymos {
namespace engine {
namespace {

std::string NormalizeCaseId(const std::string& input) {
  std::string normalized;
  normalized.reserve(input.size());
  for (size_t i = 0; i < input.size(); ++i) {
    const unsigned char ch = static_cast<unsigned char>(input[i]);
    if (std::isalnum(ch)) {
      normalized.push_back(static_cast<char>(std::tolower(ch)));
    }
  }
  return normalized;
}

}  // namespace

CaseRegistry::CaseRegistry() {
  entries_.push_back(CaseEntry{"couette", "1d", "fim", "hymos_control.sh",
                               &RunCouetteInteractive,
                               &RunCouetteBatch});
  entries_.push_back(CaseEntry{"shockstructure", "1d", "fim", "hymos_control.sh",
                               &RunShockStructureInteractive,
                               &RunShockStructureBatch});
  entries_.push_back(CaseEntry{"cavity", "2d", "fim", "hymos_control.sh",
                               &RunCavityInteractive,
                               &RunCavityBatch});
}

const CaseRegistry& CaseRegistry::Instance() {
  static const CaseRegistry registry;
  return registry;
}

const CaseEntry* CaseRegistry::Find(const std::string& case_id) const {
  const std::string normalized = NormalizeCaseId(case_id);
  for (size_t i = 0; i < entries_.size(); ++i) {
    if (NormalizeCaseId(entries_[i].id) == normalized) {
      return &entries_[i];
    }
  }
  return NULL;
}

const std::vector<CaseEntry>& CaseRegistry::List() const {
  return entries_;
}

}  // namespace engine
}  // namespace hymos
