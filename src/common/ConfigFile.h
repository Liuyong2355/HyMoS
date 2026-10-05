// HyMoS configuration reader. Independently implemented for the project format.
#ifndef HYMOS_COMMON_CONFIGFILE_H
#define HYMOS_COMMON_CONFIGFILE_H

#include <map>
#include <string>
#include "chameleon.h"

// Case-sensitive INI sections and keys. Map input uses "section/key" identifiers.
class ConfigFile {
 public:
  explicit ConfigFile(const std::string& filename);
  explicit ConfigFile(const std::map<std::string, std::string>& entries);

  // Required entries throw std::string when absent. Defaults are inserted once.
  const Chameleon& Value(const std::string& section,
                         const std::string& entry) const;
  const Chameleon& Value(const std::string& section, const std::string& entry,
                         double fallback);
  const Chameleon& Value(const std::string& section, const std::string& entry,
                         const std::string& fallback);

 private:
  std::map<std::string, Chameleon> values_;
};

#endif