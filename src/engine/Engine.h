// Usage:
// - `Engine::ListCases()` returns the cases exposed to CLI and future GUI code.
// - `RunInteractive(case_id, script)` starts the backend service for a case.
// - `RunBatch(case_id, config)` reads one config file, runs once, and exits.

#ifndef HYMOS_ENGINE_ENGINE_H
#define HYMOS_ENGINE_ENGINE_H

#include <string>
#include <vector>

#include "CaseRegistry.h"

namespace hymos {
namespace engine {

class Engine {
public:
  Engine();

  const CaseEntry* FindCase(const std::string& case_id) const;
  const std::vector<CaseEntry>& ListCases() const;

  int RunInteractive(const std::string& case_id, const std::string& control_script = "") const;
  int RunBatch(const std::string& case_id, const std::string& config_path) const;

private:
  const CaseRegistry& registry_;
};

}  // namespace engine
}  // namespace hymos

#endif
