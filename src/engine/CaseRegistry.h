// Usage:
// - Include this header when you need to discover or launch built-in cases.
// - Call `hymos::engine::CaseRegistry::Instance().List()` to enumerate cases.
// - Call `Find(case_id)` and then use `RunInteractive(...)` or `RunBatch(...)`.

#ifndef HYMOS_ENGINE_CASEREGISTRY_H
#define HYMOS_ENGINE_CASEREGISTRY_H

#include <string>
#include <vector>

namespace hymos {
namespace engine {

struct CaseEntry {
  std::string id;
  std::string dimension;
  std::string solver_family;
  std::string default_control_script;

  int (*run_interactive)(const std::string& control_script);
  int (*run_batch)(const std::string& config_path);
};

class CaseRegistry {
public:
  static const CaseRegistry& Instance();

  const CaseEntry* Find(const std::string& case_id) const;
  const std::vector<CaseEntry>& List() const;

private:
  CaseRegistry();

  std::vector<CaseEntry> entries_;
};

}  // namespace engine
}  // namespace hymos

#endif
