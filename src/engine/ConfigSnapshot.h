#ifndef HYMOS_ENGINE_CONFIGSNAPSHOT_H
#define HYMOS_ENGINE_CONFIGSNAPSHOT_H

#include <map>
#include <string>
#include <vector>

namespace hymos {
namespace engine {

// Immutable input captured before a task starts. The worker must use entries,
// never reopen source_path: callers may edit or remove the original TXT later.
struct ConfigSnapshot {
  std::string case_id;
  std::string source_path;
  std::string text;
  std::map<std::string, std::string> entries;
};

struct ConfigValidation {
  bool valid = false;
  std::string canonical_case_id;
  std::vector<std::string> errors;

  std::string Message() const;
};

// Validates a registered case and the legacy INI-like TXT grammar, then captures
// the exact text and parsed key/value map. No solver is constructed here.
// snapshot is populated only on success and is safe to pass to a background task.
ConfigValidation ValidateAndSnapshotConfig(const std::string& case_id,
                                           const std::string& config_path,
                                           ConfigSnapshot* snapshot);

// Validates legacy INI-like text and captures it exactly as a frozen snapshot.
// source_path is descriptive only: workers must never reopen it. This is the
// text counterpart of ValidateAndSnapshotConfig; both APIs use one parser.
ConfigValidation ValidateAndSnapshotConfigText(const std::string& case_id,
                                               const std::string& text,
                                               const std::string& source_path,
                                               ConfigSnapshot* snapshot);

}  // namespace engine
}  // namespace hymos

#endif
