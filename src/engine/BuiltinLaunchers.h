// Usage:
// - Declares the built-in case launchers used by the unified backend registry.
// - Keep one declaration per built-in case and implement them in dimension-
//   specific source files to avoid mixing 1D and 2D template worlds.

#ifndef HYMOS_ENGINE_BUILTINLAUNCHERS_H
#define HYMOS_ENGINE_BUILTINLAUNCHERS_H

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include <fstream>
#include <string>

namespace hymos {
namespace engine {

class RunWorkspace {
public:
  RunWorkspace(const std::string& case_id,
               const std::string& case_label,
               const std::string& config_path)
      : case_id_(case_id), case_label_(case_label), config_path_(config_path) {}

  bool Prepare() {
    old_cwd_ = CurrentDirectory();
    if (old_cwd_.empty()) {
      return false;
    }

    const std::string absolute_config = ResolveAbsolutePath(config_path_);
    config_dir_ = Dirname(absolute_config);
    config_filename_ = Basename(absolute_config);
    if (config_dir_.empty() || config_filename_.empty()) {
      return false;
    }

    const std::string runs_root = ResolveRunsRoot(config_dir_);
    if (!EnsureDirectory(runs_root)) {
      return false;
    }

    run_dir_ = runs_root + "/" + case_label_ + "_" + Timestamp();
    if (!EnsureDirectory(run_dir_)) {
      return false;
    }

    const std::string target_config = run_dir_ + "/" + config_filename_;
    if (!CopyFile(absolute_config, target_config)) {
      return false;
    }

    return chdir(run_dir_.c_str()) == 0;
  }

  ~RunWorkspace() {
    if (!old_cwd_.empty()) {
      const int restore_status = chdir(old_cwd_.c_str());
      (void)restore_status;
    }
  }

  const std::string& ConfigFilename() const { return config_filename_; }
  const std::string& RunDirectory() const { return run_dir_; }

private:
  static std::string CurrentDirectory() {
    char buffer[PATH_MAX];
    if (getcwd(buffer, sizeof(buffer)) == NULL) {
      return "";
    }
    return buffer;
  }

  static std::string ResolveAbsolutePath(const std::string& path) {
    if (!path.empty() && path[0] == '/') {
      return path;
    }
    return CurrentDirectory() + "/" + path;
  }

  static std::string Dirname(const std::string& path) {
    const std::string::size_type pos = path.find_last_of('/');
    if (pos == std::string::npos) {
      return ".";
    }
    return path.substr(0, pos);
  }

  static std::string Basename(const std::string& path) {
    const std::string::size_type pos = path.find_last_of('/');
    if (pos == std::string::npos) {
      return path;
    }
    return path.substr(pos + 1);
  }

  static std::string ResolveRunsRoot(const std::string& config_dir) {
    const char* env_root = getenv("HYMOS_RUNS_ROOT");
    if (env_root != NULL && env_root[0] != '\0') {
      return env_root;
    }
    return config_dir + "/runs";
  }

  static bool EnsureDirectory(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) == 0) {
      return S_ISDIR(st.st_mode);
    }
    return mkdir(path.c_str(), 0755) == 0;
  }

  static bool CopyFile(const std::string& src, const std::string& dst) {
    std::ifstream is(src.c_str(), std::ios::binary);
    std::ofstream os(dst.c_str(), std::ios::binary);
    os << is.rdbuf();
    return is.good() || is.eof();
  }

  static std::string Timestamp() {
    time_t now = time(NULL);
    struct tm local_tm;
    localtime_r(&now, &local_tm);
    char buffer[32];
    strftime(buffer, sizeof(buffer), "%Y%m%d_%H%M%S", &local_tm);
    return buffer;
  }

private:
  std::string case_id_;
  std::string case_label_;
  std::string config_path_;
  std::string old_cwd_;
  std::string config_dir_;
  std::string config_filename_;
  std::string run_dir_;
};

int RunCouetteInteractive(const std::string& control_script);
int RunCouetteBatch(const std::string& config_path);

int RunShockStructureInteractive(const std::string& control_script);
int RunShockStructureBatch(const std::string& config_path);

int RunCavityInteractive(const std::string& control_script);
int RunCavityBatch(const std::string& config_path);

}  // namespace engine
}  // namespace hymos

#endif
