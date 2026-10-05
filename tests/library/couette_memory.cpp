#include "Couette.h"
#include "ConfigSnapshot.h"
#include <omp.h>
#include <unistd.h>
#include <dirent.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace hymos::engine;

void Require(bool value, const char* message) {
  if (!value) throw std::runtime_error(message);
}
std::string Cwd() {
  char buf[4096];
  Require(getcwd(buf, sizeof(buf)) != NULL, "getcwd failed");
  return buf;
}
std::vector<std::string> Files() {
  DIR* dir = opendir(".");
  Require(dir != NULL, "opendir failed");
  std::vector<std::string> names;
  while (dirent* entry = readdir(dir)) names.push_back(entry->d_name);
  closedir(dir);
  std::sort(names.begin(), names.end());
  return names;
}
bool Close(double actual, double expected) {
  // Legacy files contain 12 significant digits, while memory preserves doubles.
  return std::isfinite(actual) && std::abs(actual - expected) <=
         1e-11 * std::max(1.0, std::abs(expected));
}
std::string ReadBytes(const std::string& path) {
  std::ifstream input(path.c_str(), std::ios::binary);
  Require(input.good(), "file unavailable");
  std::ostringstream data;
  data << input.rdbuf();
  Require(input.good() || input.eof(), "file read failed");
  return data.str();
}
void CompareField(const CouetteResult& result, const std::string& path) {
  std::ifstream in(path);
  Require(in.good(), "baseline field unavailable");
  std::string line;
  size_t row = 0;
  while (std::getline(in, line)) {
    std::istringstream is(line);
    double expected[11];
    if (!(is >> expected[0])) continue;
    for (int j=1; j<11; ++j) Require(bool(is >> expected[j]), "invalid baseline field");
    Require(row < result.x.size(), "field too short");
    const auto& v=result.velocity[row]; const auto& s=result.stress_moments[row];
    const auto& q=result.heat_flux[row];
    const double actual[] = {result.x[row], result.density[row], v[0], v[1], v[2],
                            result.temperature[row], s[0], s[1], s[2], q[0], q[1]};
    for (int j=0; j<11; ++j) Require(Close(actual[j], expected[j]), "field mismatch");
    ++row;
  }
  Require(row == 128 && row == result.x.size(), "unexpected field size");
}
void CompareResidual(const CouetteResult& result, const std::string& path) {
  std::ifstream in(path);
  Require(in.good(), "baseline residual unavailable");
  std::string line; size_t row=0;
  while (std::getline(in,line)) {
    std::istringstream is(line); size_t step; double value;
    if (!(is >> step >> value)) continue;
    Require(step == row && row < result.residual.size(), "residual step mismatch");
    Require(Close(result.residual[row],value), "residual mismatch"); ++row;
  }
  Require(row == result.residual.size() && row == 33, "residual history size mismatch");
}
int main(int argc, char** argv) {
  try {
    Require(argc == 5,
            "expected baseline field, residual, distribution, and input paths");
    ConfigSnapshot snapshot;
    ConfigValidation validation = ValidateAndSnapshotConfig("CoUeTtE", argv[4], &snapshot);
    Require(validation.valid, validation.Message().c_str());
    Require(snapshot.case_id == "couette" && !snapshot.text.empty(), "snapshot missing identity or text");
    Require(snapshot.entries.find("System/ORDER") != snapshot.entries.end(), "snapshot missing ORDER");
    Require(snapshot.entries["Solver/FIM"] == "3", "snapshot captured wrong FIM");
    ConfigSnapshot text_snapshot;
    ConfigValidation text_validation = ValidateAndSnapshotConfigText(
        "couette", snapshot.text, "editor://couette", &text_snapshot);
    Require(text_validation.valid, text_validation.Message().c_str());
    Require(text_snapshot.source_path == "editor://couette" &&
            text_snapshot.text == snapshot.text && text_snapshot.entries == snapshot.entries,
            "text configuration snapshot differs from file configuration snapshot");
    Require(!ValidateAndSnapshotConfigText("couette", "[System\nORDER = 7\n",
                                           "editor://bad", NULL).valid,
            "malformed text configuration accepted");
    Require(!ValidateAndSnapshotConfig("not-a-case", argv[4], NULL).valid, "unknown case accepted");
    {
      std::ofstream malformed("hymos_bad_config.txt");
      malformed << "[System\nORDER = 7\n";
    }
    ConfigValidation malformed = ValidateAndSnapshotConfig("couette", "hymos_bad_config.txt", NULL);
    Require(!malformed.valid && !malformed.errors.empty(), "malformed config accepted");
    Require(unlink("hymos_bad_config.txt") == 0, "temporary malformed config not removed");
    const auto files=Files(); const auto cwd=Cwd();
    const int previous_threads=omp_get_max_threads();
    CouetteConfig cfg;
    CouetteConfig config_from_text;
    Require(BuildCouetteConfig(snapshot, &config_from_text).valid,
            "verified TXT was rejected");
    Require(config_from_text.nx == 128 && Close(config_from_text.kn, 0.01),
            "TXT values were not captured");
    ConfigSnapshot invalid_order = snapshot;
    invalid_order.entries["System/ORDER"] = "2";
    Require(!BuildCouetteConfig(invalid_order, &config_from_text).valid,
            "ORDER below the supported moment range was accepted");
    invalid_order.entries["System/ORDER"] = "50";
    Require(!BuildCouetteConfig(invalid_order, &config_from_text).valid,
            "ORDER at the MAX_ORDER table boundary was accepted");
    ConfigSnapshot invalid_threads = snapshot;
    invalid_threads.entries["System/n_thread"] = "257";
    Require(!BuildCouetteConfig(invalid_threads, &config_from_text).valid,
            "thread count above the supported bound was accepted");
    ConfigSnapshot deferred_mesh = snapshot;
    deferred_mesh.entries["Mesh/Type"] = "1";
    Require(!BuildCouetteConfig(deferred_mesh, &config_from_text).valid,
            "deferred file-backed mesh was accepted");
    Require(config_from_text.order == 7 && config_from_text.n_thread == 1 &&
            config_from_text.mesh_type == 0 && config_from_text.solver_type == "single" &&
            config_from_text.fim == 3 && Close(config_from_text.temperature, 1.0) &&
            Close(config_from_text.ma, 1.0) && Close(config_from_text.mesh_left, -0.5) &&
            Close(config_from_text.mesh_length, 1.0), "extended TXT values were not captured");

    ConfigSnapshot expanded_snapshot = snapshot;
    expanded_snapshot.entries["System/ORDER"] = "5";
    expanded_snapshot.entries["System/n_thread"] = "2";
    expanded_snapshot.entries["System/Temperature"] = "1.1";
    expanded_snapshot.entries["Const/Ma"] = "0.8";
    expanded_snapshot.entries["Mesh/L0"] = "-1.0";
    expanded_snapshot.entries["Mesh/Len"] = "2.0";
    CouetteConfig expanded_config;
    ConfigValidation expanded_validation =
        BuildCouetteConfig(expanded_snapshot, &expanded_config);
    Require(expanded_validation.valid, expanded_validation.Message().c_str());
    Require(expanded_config.order == 5 && expanded_config.n_thread == 2 &&
            Close(expanded_config.temperature, 1.1) &&
            Close(expanded_config.ma, 0.8) &&
            Close(expanded_config.mesh_left, -1.0) &&
            Close(expanded_config.mesh_length, 2.0),
            "expanded configuration values were not preserved");
    CouetteResult expanded_result = SolveCouette(expanded_config);
    Require(expanded_result.status == CouetteStatus::Converged &&
            expanded_result.x.size() == 128 &&
            expanded_result.x.front() < -0.9 && expanded_result.x.back() > 0.9,
            "expanded ORDER/thread/physical/domain configuration failed");

    for (int fim_mode = 1; fim_mode <= 2; ++fim_mode) {
      ConfigSnapshot fim_snapshot = snapshot;
      fim_snapshot.entries["Solver/FIM"] = std::to_string(fim_mode);
      CouetteConfig fim_config;
      Require(BuildCouetteConfig(fim_snapshot, &fim_config).valid, "FIM mode was rejected");
      CouetteResult fim_result = SolveCouette(fim_config);
      Require(fim_result.status == CouetteStatus::Converged && fim_result.x.size() == 128, "FIM mode did not converge");
    }
    ConfigSnapshot nmg_snapshot = snapshot;
    nmg_snapshot.entries["Solver/SolverType"] = "NMG";
    CouetteConfig nmg_config;
    ConfigValidation nmg_validation = BuildCouetteConfig(nmg_snapshot, &nmg_config);
    Require(nmg_validation.valid, nmg_validation.Message().c_str());
    CouetteResult nmg_result = SolveCouette(nmg_config);
    Require(nmg_result.status == CouetteStatus::Converged && nmg_result.x.size() == 128, "NMG did not converge");
    nmg_snapshot.entries["Solver/FIM"] = "2";
    Require(BuildCouetteConfig(nmg_snapshot, &nmg_config).valid,
            "NMG FIM-2 configuration was rejected");
    nmg_result = SolveCouette(nmg_config);
    Require(nmg_result.status == CouetteStatus::Converged &&
            nmg_result.x.size() == 128,
            "NMG FIM-2 did not converge");
    auto a=SolveCouette(cfg);
    Require(a.status == CouetteStatus::Converged, a.message.c_str());
    CompareField(a,argv[1]); CompareResidual(a,argv[2]);
    Require(omp_get_max_threads() == previous_threads, "OpenMP setting leaked");

    auto bad=cfg; bad.kn=std::numeric_limits<double>::quiet_NaN();
    auto failed=SolveCouette(bad);
    Require(failed.status == CouetteStatus::InvalidConfig && failed.x.empty(), "NaN accepted");
    bad=cfg; bad.nx=0;
    Require(SolveCouette(bad).status == CouetteStatus::InvalidConfig, "invalid grid accepted");
    bad=cfg; bad.max_steps=0;
    Require(SolveCouette(bad).status == CouetteStatus::InvalidConfig, "unbounded solve accepted");
    bad=cfg; bad.cfl=std::numeric_limits<double>::infinity();
    Require(SolveCouette(bad).status == CouetteStatus::InvalidConfig, "infinite CFL accepted");

    auto short_cfg=cfg; short_cfg.max_steps=1;
    auto partial=SolveCouette(short_cfg);
    Require(partial.status == CouetteStatus::IterationLimit && partial.x.size()==128 &&
            partial.residual.size()==2, "iteration limit semantics incorrect");
    auto changed=cfg; changed.kn=0.02;
    auto b=SolveCouette(changed);
    Require(b.status == CouetteStatus::Converged, b.message.c_str());
    Require(a.density != b.density, "changed Kn had no effect");
    auto again=SolveCouette(cfg);
    Require(again.status == CouetteStatus::Converged, again.message.c_str());
    Require(a.x == again.x && a.density == again.density && a.velocity == again.velocity &&
            a.temperature == again.temperature && a.stress_moments == again.stress_moments &&
            a.heat_flux == again.heat_flux && a.residual == again.residual,
            "A-B-A or failure recovery changed results");
    Require(Cwd()==cwd && Files()==files, "memory solve changed cwd or files");

    ConfigValidation submitted_validation;
    auto text_task = SubmitCouette(snapshot, &submitted_validation);
    Require(submitted_validation.valid, submitted_validation.Message().c_str());
    Require(!text_task->workspace().empty() &&
            !text_task->resolved_config_path().empty(),
            "snapshot task workspace was not created");
    {
      std::ifstream frozen(text_task->resolved_config_path().c_str(),
                           std::ios::binary);
      std::ostringstream frozen_text;
      frozen_text << frozen.rdbuf();
      Require(frozen.good() || frozen.eof(),
              "frozen input.resolved.txt could not be read");
      Require(frozen_text.str() == snapshot.text,
              "frozen input.resolved.txt differs from submitted snapshot");
    }
    Require(text_task->WaitFor(std::chrono::milliseconds(60000)),
            "bounded wait did not observe completion");
    const CouetteResult text_finished = text_task->GetResult();
    Require(text_finished.status == CouetteStatus::Converged &&
            text_finished.density == a.density && text_finished.residual == a.residual,
            "TXT task differs from direct library call");
    Require(text_finished.workspace == text_task->workspace() &&
            text_finished.resolved_config_path ==
                text_task->resolved_config_path(),
            "task result lost workspace metadata");
    Require(text_finished.config.legacy_distribution_output_path ==
                text_task->workspace() + "/DisSol1th.dat",
            "task result lost legacy distribution path");
    Require(ReadBytes(text_finished.config.legacy_distribution_output_path) ==
                ReadBytes(argv[3]),
            "snapshot task distribution output differs from legacy baseline");
    const CouetteProgress final_progress = text_task->progress();
    Require(final_progress.iteration + 1 == text_finished.residual.size() &&
            Close(final_progress.residual, text_finished.residual.back()) &&
            final_progress.elapsed_seconds >= 0.0,
            "final task progress snapshot is inconsistent");
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    Require(text_task->elapsed_seconds() == final_progress.elapsed_seconds,
            "completed task elapsed time did not stay at its final iteration");

    auto task = SubmitCouette(cfg);
    const CouetteResult pending = task->GetResult();
    Require(pending.status == CouetteStatus::Pending || pending.status == CouetteStatus::Running ||
            pending.status == CouetteStatus::Converged, "invalid task initial state");
    task->Wait();
    const CouetteResult finished = task->GetResult();
    Require(finished.status == CouetteStatus::Converged, finished.message.c_str());
    Require(finished.density == a.density && finished.residual == a.residual,
            "async completed result differs from synchronous result");

    auto cancelled_task = SubmitCouette(cfg);
    cancelled_task->RequestStop();
    cancelled_task->Wait();
    const CouetteResult cancelled = cancelled_task->GetResult();
    Require(cancelled.status == CouetteStatus::Cancelled, "stop request did not cancel task");
    Require(cancelled.x.size() == 128 && !cancelled.residual.empty(),
            "cancelled task did not return a consistent partial result");

    std::ifstream sentinel(".control.1"); std::string value; std::getline(sentinel,value);
    Require(value=="999 test sentinel", "memory solve consumed control file");
    std::cout << "PASS configuration snapshot/validation, memory fields/residual/distribution, invalid parameters, iteration limit, A-B-A, "
                 "result ownership, no file effects, OpenMP restoration\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << '\n'; return 1;
  }
}
