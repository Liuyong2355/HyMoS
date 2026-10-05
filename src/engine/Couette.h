#ifndef HYMOS_ENGINE_COUETTE_H
#define HYMOS_ENGINE_COUETTE_H

#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <memory>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "ConfigSnapshot.h"

namespace hymos {
namespace engine {

// Experimental synchronous interface for the existing Couette solver.
// The legacy fixed Couette walls remain unchanged.
// Do not mix with legacy Engine calls in the same process concurrently.
struct CouetteConfig {
  // Parsed solver selection passed to the existing runtime-configured backend.
  int order = 7;
  int n_thread = 1;
  int mesh_type = 0;
  // Internal task-only destination for legacy-compatible full distribution
  // output. Empty for direct in-memory solves, which remain file-free.
  std::string legacy_distribution_output_path;
  std::string solver_type = "single";
  int fim = 3;
  std::map<std::string, std::string> text_entries;
  int nx = 128;
  double temperature = 1.0;
  double kn = 0.01;
  double pr = 1.0;
  double ma = 1.0;
  double mesh_left = -0.5;
  double mesh_length = 1.0;
  double tolerance = 1e-8;
  double cfl = 0.8;
  double cfl_euler = 0.8;
  int euler_steps = 600;
  unsigned int max_steps = 100000;
};

enum class CouetteStatus {
  Pending,
  Running,
  Converged,
  IterationLimit,
  Cancelled,
  InvalidConfig,
  Busy,
  Failed
};

struct CouetteResult {
  CouetteStatus status = CouetteStatus::Failed;
  std::string message;
  CouetteConfig config;
  // Snapshot tasks persist the exact submitted INI in their own workspace.
  std::string workspace;
  std::string resolved_config_path;
  // Cell centers on the configured domain; dimensionless variables, as in legacy output.
  std::vector<double> x, density, temperature;
  std::vector<std::array<double, 3>> velocity;
  // Legacy tauxx/tauxy/tauyy moment coefficients; not a full stress tensor.
  std::vector<std::array<double, 3>> stress_moments;
  // Legacy qx/qy components, in that order.
  std::vector<std::array<double, 2>> heat_flux;
  // PL2 residual: element 0 is initialization, followed by each outer iteration.
  std::vector<double> residual;
};

// Progress values published together at one completed outer-iteration boundary.
struct CouetteProgress {
  unsigned int iteration;
  double residual;
  double elapsed_seconds;
};

// Converts the unmodified legacy input.txt snapshot to the supported Couette
// path. Unsupported or unsafe values are reported before any task runs.
ConfigValidation BuildCouetteConfig(const ConfigSnapshot& snapshot,
                                    CouetteConfig* config);

// Owns all returned data. No config/results files or control-file polling.
// IterationLimit returns a partial field; other failures return empty fields.
CouetteResult SolveCouette(const CouetteConfig& config);
CouetteResult SolveCouette(const ConfigSnapshot& snapshot);

// A process-local task. The worker uses only libhymos; it never touches a
// Baltamatica bxArray or any other host API.
class CouetteTask {
public:
  ~CouetteTask();

  CouetteTask(const CouetteTask&) = delete;
  CouetteTask& operator=(const CouetteTask&) = delete;

  CouetteStatus status() const;
  CouetteProgress progress() const;
  unsigned int iteration() const;
  double residual() const;
  double elapsed_seconds() const;
  const std::string& workspace() const;
  const std::string& resolved_config_path() const;
  void RequestStop();
  // Returns false only when a nonnegative timeout expires while still running.
  bool WaitFor(std::chrono::milliseconds timeout);
  void Wait();
  // While pending/running, returns only status, config, and a message.
  CouetteResult GetResult() const;

private:
  explicit CouetteTask(const CouetteConfig& config,
                       const std::string& workspace = std::string(),
                       const std::string& resolved_config_path = std::string());
  void Run();

  CouetteConfig config_;
  std::string workspace_;
  std::string resolved_config_path_;
  std::atomic<CouetteStatus> status_;
  std::atomic<unsigned int> iteration_;
  std::atomic<double> residual_;
  std::atomic<long long> progress_elapsed_nanoseconds_;
  std::atomic<bool> stop_requested_;
  mutable std::mutex completion_mutex_;
  std::condition_variable completion_cv_;
  mutable std::mutex result_mutex_;
  CouetteResult result_;
  std::thread worker_;

  friend std::shared_ptr<CouetteTask> SubmitCouette(const CouetteConfig& config);
  friend std::shared_ptr<CouetteTask> SubmitCouette(const ConfigSnapshot& snapshot,
                                                    ConfigValidation* validation);
};

// Returns immediately. At most one in-memory Couette task can own the solver.
// A second submission produces a task whose final result is Busy.
std::shared_ptr<CouetteTask> SubmitCouette(const CouetteConfig& config);
// The snapshot is parsed synchronously and copied into the task before its worker
// starts. Invalid input returns an already-completed InvalidConfig task.
std::shared_ptr<CouetteTask> SubmitCouette(const ConfigSnapshot& snapshot,
                                           ConfigValidation* validation);

}  // namespace engine
}  // namespace hymos
#endif
