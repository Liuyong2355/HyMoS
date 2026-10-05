// Usage:
// - Contains launcher bindings for all built-in 2D cases.
// - Add new 2D cases here when you want them available in the unified backend.

#include "BuiltinLaunchers.h"

#include "apps/support/CaseRunner2D.h"
#include "cases/2d/CavityCase.h"

namespace hymos {
namespace engine {
namespace {

template <typename App>
int RunInteractive(const std::string& control_script) {
  App app;
  return hymos::apps::CaseRunner2D<App>::Main(app, control_script.c_str());
}

template <typename App>
int RunBatch(const std::string& case_id,
             const std::string& case_label,
             const std::string& config_path) {
  RunWorkspace workspace(case_id, case_label, config_path);
  if (!workspace.Prepare()) {
    return 1;
  }
  App app;
  app.read_config(workspace.ConfigFilename());
  app.run();
  return 0;
}

}  // namespace

int RunCavityInteractive(const std::string& control_script) {
  return RunInteractive<hymos::cases::twod::cavity::app_t>(control_script);
}

int RunCavityBatch(const std::string& config_path) {
  return RunBatch<hymos::cases::twod::cavity::app_t>("cavity", "Cavity", config_path);
}

}  // namespace engine
}  // namespace hymos
