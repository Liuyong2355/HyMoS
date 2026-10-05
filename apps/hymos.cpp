// Usage:
// - `hymos list-cases`
// - `hymos serve <case-id> [control-script]`
// - `hymos run <case-id> <config-file>`
// This executable is the unified backend entry for CLI, scripts, and future GUI.

#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "Engine.h"

namespace {

void PrintUsage(std::ostream& os) {
  os << "Usage:\n";
  os << "  hymos list-cases\n";
  os << "  hymos serve <case-id> [control-script]\n";
  os << "  hymos run <case-id> <config-file>\n";
}

int ListCases() {
  hymos::engine::Engine engine;
  const std::vector<hymos::engine::CaseEntry>& entries = engine.ListCases();
  for (size_t i = 0; i < entries.size(); ++i) {
    std::cout << entries[i].id << "\t"
              << entries[i].dimension << "\t"
              << entries[i].solver_family << "\t"
              << entries[i].default_control_script << std::endl;
  }
  return 0;
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    if (argc < 2) {
      PrintUsage(std::cerr);
      return 1;
    }

    const std::string command(argv[1]);
    hymos::engine::Engine engine;

    if (command == "list-cases") {
      return ListCases();
    }

    if (command == "serve") {
      if (argc < 3) {
        PrintUsage(std::cerr);
        return 1;
      }
      const std::string case_id(argv[2]);
      const std::string script_name = (argc >= 4) ? argv[3] : "";
      return engine.RunInteractive(case_id, script_name);
    }

    if (command == "run") {
      if (argc < 4) {
        PrintUsage(std::cerr);
        return 1;
      }
      return engine.RunBatch(argv[2], argv[3]);
    }

    PrintUsage(std::cerr);
    return 1;
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  }
}
