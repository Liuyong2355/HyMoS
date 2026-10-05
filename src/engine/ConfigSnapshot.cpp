#include "ConfigSnapshot.h"

#include <fstream>
#include <sstream>

#include "CaseRegistry.h"

namespace hymos {
namespace engine {
namespace {

std::string Trim(const std::string& input) {
  const std::string whitespace(" \t\r\n");
  const std::string::size_type first = input.find_first_not_of(whitespace);
  if (first == std::string::npos) return "";
  const std::string::size_type last = input.find_last_not_of(whitespace);
  return input.substr(first, last - first + 1);
}

void AddError(ConfigValidation* validation, size_t line, const std::string& message) {
  std::ostringstream os;
  if (line != 0) os << "line " << line << ": ";
  os << message;
  validation->errors.push_back(os.str());
}

}  // namespace

std::string ConfigValidation::Message() const {
  if (errors.empty()) return "OK";
  std::ostringstream os;
  for (size_t i = 0; i < errors.size(); ++i) {
    if (i != 0) os << "; ";
    os << errors[i];
  }
  return os.str();
}

ConfigValidation ValidateAndSnapshotConfig(const std::string& case_id,
                                           const std::string& config_path,
                                           ConfigSnapshot* snapshot) {
  if (config_path.empty()) {
    ConfigValidation validation;
    validation.errors.push_back("configuration path is empty");
    return validation;
  }

  std::ifstream input(config_path.c_str(), std::ios::binary);
  if (!input.good()) {
    ConfigValidation validation;
    validation.errors.push_back("cannot open configuration file: " + config_path);
    return validation;
  }
  std::ostringstream buffer;
  buffer << input.rdbuf();
  if (input.bad()) {
    ConfigValidation validation;
    validation.errors.push_back("cannot read configuration file: " + config_path);
    return validation;
  }
  return ValidateAndSnapshotConfigText(case_id, buffer.str(), config_path, snapshot);
}

ConfigValidation ValidateAndSnapshotConfigText(const std::string& case_id,
                                               const std::string& text,
                                               const std::string& source_path,
                                               ConfigSnapshot* snapshot) {
  ConfigValidation validation;
  const CaseEntry* entry = CaseRegistry::Instance().Find(case_id);
  if (entry == NULL) {
    validation.errors.push_back("unknown case: " + case_id);
    return validation;
  }
  validation.canonical_case_id = entry->id;
  ConfigSnapshot candidate;
  candidate.case_id = entry->id;
  candidate.source_path = source_path;
  candidate.text = text;

  std::istringstream lines(candidate.text);
  std::string raw;
  std::string section;
  size_t line_number = 0;
  while (std::getline(lines, raw)) {
    ++line_number;
    const std::string line = Trim(raw);
    if (line.empty() || line[0] == '#' || line[0] == ';') continue;
    if (line[0] == '[') {
      const std::string::size_type close = line.find(']');
      if (close == std::string::npos || close != line.size() - 1) {
        AddError(&validation, line_number, "section header must be [name]");
        continue;
      }
      section = Trim(line.substr(1, close - 1));
      if (section.empty()) AddError(&validation, line_number, "section name is empty");
      continue;
    }
    const std::string::size_type equals = line.find('=');
    if (equals == std::string::npos) {
      AddError(&validation, line_number, "expected key=value");
      continue;
    }
    if (section.empty()) {
      AddError(&validation, line_number, "key/value appears before a section");
      continue;
    }
    const std::string key = Trim(line.substr(0, equals));
    const std::string value = Trim(line.substr(equals + 1));
    if (key.empty()) {
      AddError(&validation, line_number, "key is empty");
      continue;
    }
    candidate.entries[section + "/" + key] = value;
  }

  if (candidate.entries.empty() && validation.errors.empty()) {
    validation.errors.push_back("configuration contains no key/value entries");
  }
  if (!validation.errors.empty()) return validation;

  validation.valid = true;
  if (snapshot != NULL) *snapshot = candidate;
  return validation;
}

}  // namespace engine
}  // namespace hymos
