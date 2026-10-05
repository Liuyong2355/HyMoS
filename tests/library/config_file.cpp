#include "ConfigFile.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

template <typename Operation>
void expect_string_error(Operation operation, const std::string& expected) {
  try {
    operation();
  } catch (const std::string& message) {
    require(message == expected, "Unexpected configuration error: " + message);
    return;
  }
  throw std::runtime_error("Missing configuration error: " + expected);
}

struct TemporaryInput {
  std::string directory;
  std::string filename;

  TemporaryInput() {
    char pattern[] = "/tmp/hymos-config-contract-XXXXXX";
    const char* created = mkdtemp(pattern);
    if (!created) throw std::runtime_error("Cannot create temporary directory");
    directory = created;
    filename = directory + "/input.txt";
  }

  void write(const std::string& text) const {
    std::ofstream output(filename.c_str(), std::ios::binary | std::ios::trunc);
    output << text;
    if (!output.good()) throw std::runtime_error("Cannot write test input");
  }

  ~TemporaryInput() {
    std::remove(filename.c_str());
    rmdir(directory.c_str());
  }
};

void check_map_and_defaults() {
  const std::map<std::string, std::string> entries = {
      {"System/ORDER", "7"}, {"Const/Kn", "0.01"}, {"Text/value", "a=b"}};
  ConfigFile config(entries);
  const ConfigFile& required = config;
  require(static_cast<int>(required.Value("System", "ORDER")) == 7,
          "Map integer value");
  require(std::abs(static_cast<double>(required.Value("Const", "Kn")) - 0.01)
              < 1e-14, "Map floating-point value");
  require(static_cast<std::string>(required.Value("Text", "value")) == "a=b",
          "Map string value");
  expect_string_error([&] { required.Value("System", "order"); },
                      "System/order does not exist");
  require(static_cast<int>(config.Value("System", "ORDER", 99.0)) == 7,
          "Default must preserve existing numeric value");
  const Chameleon& fallback = config.Value("System", "n_thread", 4.0);
  config.Value("Text", "new", std::string("default"));
  require(static_cast<int>(fallback) == 4, "Default reference stays usable");
  require(static_cast<int>(required.Value("System", "n_thread")) == 4,
          "Numeric default persists");
  require(static_cast<std::string>(config.Value("Text", "new", std::string("other")))
              == "default", "String default persists");
}

void check_file_format(TemporaryInput& temporary) {
  temporary.write("# comment\r\n  ; indented comment\r\n\r\n"
                  " [ System ] \r\n ORDER = 3\r\nORDER=7\r\n"
                  " value = left=right=tail  \r\nempty=\r\n"
                  " literal = data # keep ; keep\r\n"
                  "[Other]\r\nORDER=9\r\n");
  const ConfigFile config(temporary.filename);
  require(static_cast<int>(config.Value("System", "ORDER")) == 7,
          "Whitespace, CRLF, sections and duplicate keys");
  require(static_cast<int>(config.Value("Other", "ORDER")) == 9,
          "Keys remain section-specific");
  require(static_cast<std::string>(config.Value("System", "value")) ==
              "left=right=tail", "Embedded equals signs");
  require(static_cast<std::string>(config.Value("System", "empty")).empty(),
          "Empty values");
  require(static_cast<std::string>(config.Value("System", "literal")) ==
              "data # keep ; keep", "Values retain comment characters");
  expect_string_error([&] { config.Value("Missing", "key"); },
                      "Missing/key does not exist");

  for (const std::string& malformed : {
           "missing_equals", "=value", "[unclosed", "[]", "[S]extra"}) {
    temporary.write("# first line\n" + malformed + "\n");
    expect_string_error([&] { ConfigFile invalid(temporary.filename); },
                        temporary.filename + ":2: invalid configuration line");
  }
  const std::string absent = temporary.directory + "/absent.txt";
  expect_string_error([&] { ConfigFile missing(absent); }, absent + " does not exist");
}

}  // namespace

int main() {
  try {
    TemporaryInput temporary;
    check_map_and_defaults();
    check_file_format(temporary);
    std::cout << "PASS configuration input contract\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
  } catch (const std::string& error) {
    std::cerr << error << '\n';
  }
  return 1;
}