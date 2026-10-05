// HyMoS configuration reader. Independently implemented for the project format.
#include "ConfigFile.h"

#include <cctype>
#include <fstream>

namespace {

std::string strip_space(const std::string& text) {
  std::size_t first = 0;
  std::size_t last = text.size();
  while (first < last && std::isspace(static_cast<unsigned char>(text[first]))) {
    ++first;
  }
  while (last > first && std::isspace(static_cast<unsigned char>(text[last - 1]))) {
    --last;
  }
  return text.substr(first, last - first);
}

std::string entry_id(const std::string& section, const std::string& entry) {
  return section + "/" + entry;
}

}  // namespace

ConfigFile::ConfigFile(const std::string& filename) {
  std::ifstream input(filename.c_str());
  if (!input.is_open()) {
    throw filename + " does not exist";
  }

  std::string section;
  std::string record;
  std::size_t line_number = 0;
  while (std::getline(input, record)) {
    ++line_number;
    const std::string text = strip_space(record);
    if (text.empty() || text.front() == '#' || text.front() == ';') {
      continue;
    }

    const auto invalid_line = [&]() -> std::string {
      return filename + ":" + std::to_string(line_number) +
             ": invalid configuration line";
    };

    if (text.front() == '[') {
      if (text.size() < 3 || text.back() != ']' ||
          text.find(']') != text.size() - 1) {
        throw invalid_line();
      }
      section = strip_space(text.substr(1, text.size() - 2));
      if (section.empty() || section.find('[') != std::string::npos) {
        throw invalid_line();
      }
      continue;
    }

    const std::size_t separator = text.find('=');
    if (separator == std::string::npos) {
      throw invalid_line();
    }
    const std::string key = strip_space(text.substr(0, separator));
    if (key.empty()) {
      throw invalid_line();
    }
    values_[entry_id(section, key)] =
        strip_space(text.substr(separator + 1));
  }
  if (input.bad()) {
    throw filename + ": unable to read configuration";
  }
}

ConfigFile::ConfigFile(const std::map<std::string, std::string>& entries) {
  for (const auto& item : entries) {
    values_.emplace(item.first, Chameleon(item.second));
  }
}

const Chameleon& ConfigFile::Value(const std::string& section,
                                 const std::string& entry) const {
  const std::string id = entry_id(section, entry);
  const auto found = values_.find(id);
  if (found == values_.end()) {
    throw id + " does not exist";
  }
  return found->second;
}

const Chameleon& ConfigFile::Value(const std::string& section,
                                 const std::string& entry, double fallback) {
  return values_.emplace(entry_id(section, entry), Chameleon(fallback))
      .first->second;
}

const Chameleon& ConfigFile::Value(const std::string& section,
                                 const std::string& entry,
                                 const std::string& fallback) {
  return values_.emplace(entry_id(section, entry), Chameleon(fallback))
      .first->second;
}