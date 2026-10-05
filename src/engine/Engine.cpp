// Usage:
// - This is the backend entry used by the new `hymos` executable.
// - Frontends should call this layer instead of directly instantiating case apps.

#include "Engine.h"

#include <stdexcept>

namespace hymos {
namespace engine {
namespace {

std::string ResolveScriptName(const CaseEntry& entry, const std::string& requested_script) {
  if (!requested_script.empty()) {
    return requested_script;
  }
  return entry.default_control_script;
}

}  // namespace

Engine::Engine() : registry_(CaseRegistry::Instance()) {}

const CaseEntry* Engine::FindCase(const std::string& case_id) const {
  return registry_.Find(case_id);
}

const std::vector<CaseEntry>& Engine::ListCases() const {
  return registry_.List();
}

int Engine::RunInteractive(const std::string& case_id, const std::string& control_script) const {
  const CaseEntry* entry = FindCase(case_id);
  if (entry == NULL) {
    throw std::runtime_error("unknown case: " + case_id);
  }
  return entry->run_interactive(ResolveScriptName(*entry, control_script));
}

int Engine::RunBatch(const std::string& case_id, const std::string& config_path) const {
  const CaseEntry* entry = FindCase(case_id);
  if (entry == NULL) {
    throw std::runtime_error("unknown case: " + case_id);
  }
  return entry->run_batch(config_path);
}

}  // namespace engine
}  // namespace hymos
