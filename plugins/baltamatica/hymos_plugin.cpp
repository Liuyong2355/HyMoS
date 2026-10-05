#include <bex/bex.hpp>
#include <dlfcn.h>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "ConfigSnapshot.h"
#include "Couette.h"

namespace {
using namespace hymos::engine;

// The 2025 Linux host owns this pointer through bxCreateCStruct. Do not use
// the C++ external-object helpers: their registration ABI is not compatible.
struct Task {
  explicit Task(const std::shared_ptr<CouetteTask>& task) : value(task) {}
  std::shared_ptr<CouetteTask> value;
};

int g_task_struct_id = 0;
bool g_task_struct_registered = false;
std::mutex g_tasks_mutex;
std::mutex g_export_mutex;
std::vector<std::weak_ptr<CouetteTask> > g_tasks;

void TrackTask(const std::shared_ptr<CouetteTask>& task) {
  std::lock_guard<std::mutex> lock(g_tasks_mutex);
  std::vector<std::weak_ptr<CouetteTask> >::iterator out = g_tasks.begin();
  for (std::vector<std::weak_ptr<CouetteTask> >::iterator it = g_tasks.begin();
       it != g_tasks.end(); ++it) {
    if (!it->expired()) *out++ = *it;
  }
  g_tasks.erase(out, g_tasks.end());
  g_tasks.push_back(task);
}

// These callbacks intentionally use no Baltamatica SDK API. A copied host
// value shares the asynchronous task; final release stops and joins it.
void* CopyTask(const void* raw) {
  const Task* task = static_cast<const Task*>(raw);
  return task ? new Task(task->value) : NULL;
}
void DeleteTask(void* raw) {
  Task* task = static_cast<Task*>(raw);
  if (!task) return;
  if (task->value && task->value.use_count() == 1) {
    task->value->RequestStop();
    task->value->Wait();
  }
  delete task;
}

const char* State(CouetteStatus state) {
  switch (state) {
    case CouetteStatus::Pending: return "pending";
    case CouetteStatus::Running: return "running";
    case CouetteStatus::Converged: return "converged";
    case CouetteStatus::IterationLimit: return "iteration_limit";
    case CouetteStatus::Cancelled: return "stopped";
    case CouetteStatus::InvalidConfig: return "invalid_config";
    case CouetteStatus::Busy: return "busy";
    default: return "failed";
  }
}
bool Done(CouetteStatus state) {
  return state != CouetteStatus::Pending && state != CouetteStatus::Running;
}
void Need(bool ok, const char* why) { if (!ok) throw std::runtime_error(why); }
std::string Text(const bxArray* input, const char* name) {
  Need(input && (bxIsString(input) || bxIsChar(input)), name);
  return baltam::bxAsString(input);
}
Task* GetTask(const bxArray* input) {
  Need(g_task_struct_registered && bxIsExternID(input, g_task_struct_id),
       "expected a HyMoS task");
  Task* task = static_cast<Task*>(bxGetCStruct(g_task_struct_id, input));
  if (!task || !task->value) throw std::runtime_error("invalid HyMoS task");
  return task;
}
bxArray* Status(const std::shared_ptr<CouetteTask>& task) {
  const char* fields[]={"state","iteration","residual","elapsed_seconds",
                       "message","workspace","resolved_config","tolerance"};
  bxArray* out=bxCreateStructMatrix(1,1,8,fields);
  CouetteResult result=task->GetResult();
  const CouetteProgress progress=task->progress();
  bxSetField(out,0,"state",bxCreateStringScalar(State(task->status())));
  bxSetField(out,0,"iteration",bxCreateDoubleScalar(progress.iteration));
  bxSetField(out,0,"residual",bxCreateDoubleScalar(progress.residual));
  bxSetField(out,0,"elapsed_seconds",bxCreateDoubleScalar(progress.elapsed_seconds));
  bxSetField(out,0,"tolerance",bxCreateDoubleScalar(result.config.tolerance));
  bxSetField(out,0,"message",bxCreateStringScalar(result.message.c_str()));
  bxSetField(out,0,"workspace",bxCreateStringScalar(task->workspace().c_str()));
  bxSetField(out,0,"resolved_config",
             bxCreateStringScalar(task->resolved_config_path().c_str()));
  return out;
}
bxArray* Column(const std::vector<double>& values) {
  bxArray* out=bxCreateDoubleMatrix((baSize)values.size(),1,bxREAL);
  double* p=bxGetDoublesRW(out);
  for(size_t i=0;i<values.size();++i) p[i]=values[i];
  return out;
}
template <size_t N> bxArray* Rows(const std::vector<std::array<double,N> >& values) {
  bxArray* out=bxCreateDoubleMatrix((baSize)values.size(),(baSize)N,bxREAL);
  double* p=bxGetDoublesRW(out);
  for(size_t i=0;i<values.size();++i) for(size_t j=0;j<N;++j)
    p[i+values.size()*j]=values[i][j];
  return out;
}
bxArray* Result(const CouetteResult& r, const CouetteProgress& progress) {
  const char* fields[]={"x","density","velocity","temperature","stress_moments","heat_flux","residual","metadata"};
  bxArray* out=bxCreateStructMatrix(1,1,8,fields);
  bxSetField(out,0,"x",Column(r.x)); bxSetField(out,0,"density",Column(r.density));
  bxSetField(out,0,"velocity",Rows(r.velocity)); bxSetField(out,0,"temperature",Column(r.temperature));
  bxSetField(out,0,"stress_moments",Rows(r.stress_moments)); bxSetField(out,0,"heat_flux",Rows(r.heat_flux));
  bxSetField(out,0,"residual",Column(r.residual));
  const char* meta_fields[]={
      "case", "dimension", "state", "message", "iteration",
      "elapsed_seconds", "order", "n_thread", "mesh_type", "nx",
      "mesh_left", "mesh_length", "temperature", "kn", "pr", "ma",
      "tolerance", "cfl", "cfl_euler", "euler_steps", "max_steps",
      "solver_type", "fim", "velocity_components", "stress_definition",
      "heat_flux_components", "field_columns", "nondimensional_definition",
      "workspace", "resolved_config", "distribution_file"};
  bxArray* meta=bxCreateStructMatrix(1,1,31,meta_fields);
  const CouetteConfig& config = r.config;
  bxSetField(meta,0,"case",bxCreateStringScalar("couette"));
  bxSetField(meta,0,"dimension",bxCreateDoubleScalar(1));
  bxSetField(meta,0,"state",bxCreateStringScalar(State(r.status)));
  bxSetField(meta,0,"message",bxCreateStringScalar(r.message.c_str()));
  bxSetField(meta,0,"iteration",bxCreateDoubleScalar(progress.iteration));
  bxSetField(meta,0,"elapsed_seconds",bxCreateDoubleScalar(progress.elapsed_seconds));
  bxSetField(meta,0,"order",bxCreateDoubleScalar(config.order));
  bxSetField(meta,0,"n_thread",bxCreateDoubleScalar(config.n_thread));
  bxSetField(meta,0,"mesh_type",bxCreateDoubleScalar(config.mesh_type));
  bxSetField(meta,0,"nx",bxCreateDoubleScalar(config.nx));
  bxSetField(meta,0,"mesh_left",bxCreateDoubleScalar(config.mesh_left));
  bxSetField(meta,0,"mesh_length",bxCreateDoubleScalar(config.mesh_length));
  bxSetField(meta,0,"temperature",bxCreateDoubleScalar(config.temperature));
  bxSetField(meta,0,"kn",bxCreateDoubleScalar(config.kn));
  bxSetField(meta,0,"pr",bxCreateDoubleScalar(config.pr));
  bxSetField(meta,0,"ma",bxCreateDoubleScalar(config.ma));
  bxSetField(meta,0,"tolerance",bxCreateDoubleScalar(config.tolerance));
  bxSetField(meta,0,"cfl",bxCreateDoubleScalar(config.cfl));
  bxSetField(meta,0,"cfl_euler",bxCreateDoubleScalar(config.cfl_euler));
  bxSetField(meta,0,"euler_steps",bxCreateDoubleScalar(config.euler_steps));
  bxSetField(meta,0,"max_steps",bxCreateDoubleScalar(config.max_steps));
  bxSetField(meta,0,"solver_type",bxCreateStringScalar(config.solver_type.c_str()));
  bxSetField(meta,0,"fim",bxCreateDoubleScalar(config.fim));
  bxSetField(meta,0,"velocity_components",bxCreateStringScalar("u1,u2,u3"));
  bxSetField(meta,0,"stress_definition",bxCreateStringScalar("legacy tauxx,tauxy,tauyy moment coefficients"));
  bxSetField(meta,0,"heat_flux_components",bxCreateStringScalar("qx,qy"));
  bxSetField(meta,0,"field_columns",bxCreateStringScalar("x,density,u1,u2,u3,temperature,tauxx,tauxy,tauyy,qx,qy"));
  bxSetField(meta,0,"nondimensional_definition",bxCreateStringScalar("legacy Couette nondimensional variables; see resolved_config and case documentation"));
  bxSetField(meta,0,"workspace",bxCreateStringScalar(r.workspace.c_str()));
  bxSetField(meta,0,"resolved_config",
             bxCreateStringScalar(r.resolved_config_path.c_str()));
  bxSetField(meta,0,"distribution_file",
             bxCreateStringScalar(r.config.legacy_distribution_output_path.c_str()));
  bxSetField(out,0,"metadata",meta); return out;
}

void CheckResultShape(const CouetteResult& result) {
  const size_t count = result.x.size();
  Need(count != 0, "task has no field data to export");
  Need(result.density.size() == count && result.velocity.size() == count &&
       result.temperature.size() == count &&
       result.stress_moments.size() == count &&
       result.heat_flux.size() == count,
       "task result arrays have inconsistent lengths");
}

void WriteAtomic(const std::string& path,
                 const std::function<void(std::ostream&)>& writer) {
  const std::string temporary = path + ".tmp";
  std::ofstream output(temporary.c_str(),
                       std::ios::out | std::ios::trunc | std::ios::binary);
  if (!output) throw std::runtime_error("cannot create export file: " + temporary);
  try {
    writer(output);
    output.close();
    if (!output) throw std::runtime_error("cannot complete export file: " + temporary);
  } catch (...) {
    output.close();
    std::remove(temporary.c_str());
    throw;
  }
  if (std::rename(temporary.c_str(), path.c_str()) != 0) {
    const std::string reason = std::strerror(errno);
    std::remove(temporary.c_str());
    throw std::runtime_error("cannot publish export file " + path + ": " + reason);
  }
}

std::string OneLine(std::string value) {
  for (size_t i = 0; i < value.size(); ++i) {
    if (value[i] == '\n' || value[i] == '\r') value[i] = ' ';
  }
  return value;
}

bxArray* ExportResultFiles(const std::shared_ptr<CouetteTask>& task) {
  const CouetteStatus state = task->status();
  Need(Done(state), "task is still running; use hymos_status or hymos_wait");
  CouetteResult result = task->GetResult();
  if (state == CouetteStatus::InvalidConfig || state == CouetteStatus::Busy ||
      state == CouetteStatus::Failed) {
    throw std::runtime_error(result.message);
  }
  CheckResultShape(result);
  Need(!result.workspace.empty(), "task has no workspace");

  const CouetteProgress progress = task->progress();
  const std::string fields_path = result.workspace + "/couette_fields.csv";
  const std::string residual_path = result.workspace + "/residual.csv";
  const std::string metadata_path = result.workspace + "/result_metadata.txt";
  const std::string distribution_path = result.workspace + "/DisSol1th.dat";
  std::lock_guard<std::mutex> lock(g_export_mutex);
  std::ifstream distribution(distribution_path.c_str(), std::ios::binary);
  Need(distribution.good(), "legacy distribution output is unavailable");

  WriteAtomic(fields_path, [&result](std::ostream& output) {
    output << "x,density,u1,u2,u3,temperature,tauxx,tauxy,tauyy,qx,qy\n";
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (size_t i = 0; i < result.x.size(); ++i) {
      output << result.x[i] << ',' << result.density[i] << ','
             << result.velocity[i][0] << ',' << result.velocity[i][1] << ','
             << result.velocity[i][2] << ',' << result.temperature[i] << ','
             << result.stress_moments[i][0] << ','
             << result.stress_moments[i][1] << ','
             << result.stress_moments[i][2] << ','
             << result.heat_flux[i][0] << ',' << result.heat_flux[i][1] << '\n';
    }
  });
  WriteAtomic(residual_path, [&result](std::ostream& output) {
    output << "iteration,residual\n";
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (size_t i = 0; i < result.residual.size(); ++i)
      output << i << ',' << result.residual[i] << '\n';
  });
  WriteAtomic(metadata_path, [&result, &progress, &distribution_path](std::ostream& output) {
    const CouetteConfig& config = result.config;
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << "case=couette\n"
           << "state=" << State(result.status) << '\n'
           << "message=" << OneLine(result.message) << '\n'
           << "iteration=" << progress.iteration << '\n'
           << "elapsed_seconds=" << progress.elapsed_seconds << '\n'
           << "order=" << config.order << '\n'
           << "n_thread=" << config.n_thread << '\n'
           << "mesh_type=" << config.mesh_type << '\n'
           << "nx=" << config.nx << '\n'
           << "mesh_left=" << config.mesh_left << '\n'
           << "mesh_length=" << config.mesh_length << '\n'
           << "temperature=" << config.temperature << '\n'
           << "kn=" << config.kn << '\n'
           << "pr=" << config.pr << '\n'
           << "ma=" << config.ma << '\n'
           << "tolerance=" << config.tolerance << '\n'
           << "cfl=" << config.cfl << '\n'
           << "cfl_euler=" << config.cfl_euler << '\n'
           << "euler_steps=" << config.euler_steps << '\n'
           << "max_steps=" << config.max_steps << '\n'
           << "solver_type=" << config.solver_type << '\n'
           << "fim=" << config.fim << '\n'
           << "distribution_file=" << distribution_path << '\n'
           << "distribution_format=legacy DisSol1th.dat internal storage order\n"
           << "field_columns=x,density,u1,u2,u3,temperature,tauxx,tauxy,tauyy,qx,qy\n"
           << "stress_definition=legacy tauxx,tauxy,tauyy moment coefficients\n"
           << "heat_flux_components=qx,qy\n"
           << "resolved_config=" << result.resolved_config_path << '\n';
  });

  const char* fields[] = {"workspace", "fields", "residual", "metadata", "distribution"};
  bxArray* summary = bxCreateStructMatrix(1, 1, 5, fields);
  bxSetField(summary, 0, "workspace",
             bxCreateStringScalar(result.workspace.c_str()));
  bxSetField(summary, 0, "fields", bxCreateStringScalar(fields_path.c_str()));
  bxSetField(summary, 0, "residual", bxCreateStringScalar(residual_path.c_str()));
  bxSetField(summary, 0, "metadata", bxCreateStringScalar(metadata_path.c_str()));
  bxSetField(summary, 0, "distribution",
             bxCreateStringScalar(distribution_path.c_str()));
  return summary;
}

bxArray* CreateTaskObject(const ConfigSnapshot& snapshot) {
  ConfigValidation validation;
  std::shared_ptr<CouetteTask> job=SubmitCouette(snapshot,&validation);
  if(!validation.valid) throw std::runtime_error(validation.Message());
  Need(g_task_struct_registered, "HyMoS task type is not initialized");
  Task* task = new Task(job);
  bxArray* object = bxCreateCStruct(g_task_struct_id, task);
  if (!object) {
    delete task;
    throw std::runtime_error("could not create HyMoS task object");
  }
  TrackTask(job);
  return object;
}

void Fail(const std::exception& e) { bxErrMsgTxt(e.what()); }


std::string PluginDirectory() {
  Dl_info info;
  if (dladdr(reinterpret_cast<void*>(&PluginDirectory), &info) == 0 ||
      info.dli_fname == NULL) {
    throw std::runtime_error("cannot locate the loaded HyMoS plugin");
  }
  const std::string path(info.dli_fname);
  const std::string::size_type slash = path.find_last_of('/');
  if (slash == std::string::npos) {
    throw std::runtime_error("invalid HyMoS plugin path: " + path);
  }
  return path.substr(0, slash);
}

void CallHostNoOutput(const char* function_name,
                      int nrhs,
                      const bxArray* inputs[]) {
  bxArray* ignored[1] = {NULL};
  const int status = bxCallBaltamatica(0, ignored, nrhs, inputs, function_name);
  if (status != 0) {
    throw std::runtime_error(std::string("Baltamatica call failed: ") + function_name);
  }
}

// The convenience entry executes only host-side script code. Submission and
// monitoring continue to use the existing public functions and task CStruct.
BALTAM_PLUGIN_FCN(HymosRun) {
  try {
    Need(nrhs == 0 && nlhs == 0, "hymos_run expects no inputs and no outputs");
    bxArray* ui_path = bxCreateStringScalar((PluginDirectory() + "/ui").c_str());
    const bxArray* inputs[1] = {ui_path};
    const int added = bxCallBaltamatica(0, NULL, 1, inputs, "addpath");
    bxDestroyArray(ui_path);
    Need(added == 0, "cannot add the HyMoS helper path");
    // Preserve the host's script error / interrupt instead of replacing it
    // with a generic C++ error after the foreground monitor is interrupted.
    bxCallBaltamatica(0, NULL, 0, NULL, "hymos_run_impl");
  } catch (const std::exception& e) {
    Fail(e);
  }
}

// Host-side editor file adapter. Resolve links before protecting installed
// packages; publish a complete file with rename rather than truncating it.
std::string ParameterFilePath(const std::string& requested) {
  Need(!requested.empty() && requested.find('\0') == std::string::npos,
       "parameter file path must be nonempty text without NUL");
  char* resolved = realpath(requested.c_str(), NULL);
  std::string path;
  if (resolved) {
    path = resolved;
    std::free(resolved);
    struct stat st;
    Need(stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode),
         "parameter destination must be a regular file");
  } else {
    Need(errno == ENOENT, "cannot resolve parameter file path");
    const size_t slash = requested.find_last_of('/');
    const std::string parent = slash == std::string::npos ? "."
        : (slash == 0 ? "/" : requested.substr(0, slash));
    const std::string name = slash == std::string::npos ? requested
        : requested.substr(slash + 1);
    Need(!name.empty() && name != "." && name != "..", "invalid parameter filename");
    resolved = realpath(parent.c_str(), NULL);
    Need(resolved != NULL, "parameter destination directory does not exist");
    path = std::string(resolved) + "/" + name;
    std::free(resolved);
    struct stat st;
    // Do not replace dangling symlinks or other unresolved existing objects.
    Need(lstat(path.c_str(), &st) != 0 && errno == ENOENT,
         "cannot save through an unresolved parameter file link");
  }
  // Keep the loaded package read-only wherever the user installed it; do not
  // assume the standard /opt/Baltamatica/plugins location.
  const std::string protected_root = PluginDirectory();
  Need(path != protected_root && path.compare(0, protected_root.size() + 1,
       protected_root + "/") != 0, "installed plugin files are read-only here; save to a user file");
  return path;
}

void SaveParameterText(const std::string& path, const std::string& text) {
  struct stat old;
  const bool existed = stat(path.c_str(), &old) == 0;
  if (existed) Need(access(path.c_str(), W_OK) == 0, "parameter file is not writable");
  std::string pattern = path + ".hymos-XXXXXX";
  std::vector<char> name(pattern.begin(), pattern.end());
  name.push_back('\0');
  int fd = mkstemp(name.data());
  Need(fd >= 0, "cannot create a temporary parameter file beside the destination");
  try {
    if (existed) Need(fchmod(fd, old.st_mode & 0777) == 0, "cannot preserve parameter permissions");
    size_t offset = 0;
    while (offset < text.size()) {
      const ssize_t count = write(fd, text.data() + offset, text.size() - offset);
      if (count < 0 && errno == EINTR) continue;
      Need(count > 0, "cannot write complete parameter text");
      offset += static_cast<size_t>(count);
    }
    Need(fsync(fd) == 0, "cannot flush parameter file");
    const int closed = close(fd);
    fd = -1;
    Need(closed == 0, "cannot close parameter file");
    Need(std::rename(name.data(), path.c_str()) == 0, "cannot replace parameter file");
  } catch (...) {
    if (fd >= 0) close(fd);
    std::remove(name.data());
    throw;
  }
}

BALTAM_PLUGIN_FCN(HymosUIFile) {
  try {
    Need((nrhs == 1 || nrhs == 2) && nlhs == 1,
         "hymos_ui_file expects path, optional text, and one path output");
    const std::string path = ParameterFilePath(Text(prhs[0], "file path must be text"));
    if (nrhs == 2) {
      const std::string text = Text(prhs[1], "configuration must be text");
      ConfigSnapshot snapshot;
      ConfigValidation check = ValidateAndSnapshotConfigText("couette", text, path, &snapshot);
      if (!check.valid) throw std::runtime_error(check.Message());
      CouetteConfig config;
      check = BuildCouetteConfig(snapshot, &config);
      if (!check.valid) throw std::runtime_error(check.Message());
      SaveParameterText(path, text);
    }
    plhs[0] = bxCreateStringScalar(path.c_str());
  } catch (const std::exception& e) { Fail(e); }
}

// Lossless form/text bridge. The original snapshot parser remains authoritative.
// Only the last occurrence of an edited key is changed (the parser's semantics).
std::string FormTrim(const std::string& s) {
  const size_t first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  return s.substr(first, s.find_last_not_of(" \t\r\n") - first + 1);
}
std::string PatchFormText(const std::string& text, const std::string& key,
                          const std::string& value) {
  const size_t split = key.find('/');
  Need(split != std::string::npos && split > 0 && split + 1 < key.size() &&
       key.find_first_of("\r\n=[]") == std::string::npos,
       "form key must be Section/Key");
  const std::string target_section = key.substr(0, split);
  const std::string target_name = key.substr(split + 1);
  std::string section;
  size_t value_begin = std::string::npos, value_end = 0;
  size_t insertion = std::string::npos;
  bool found_section = false;
  for (size_t pos = 0; pos < text.size();) {
    const size_t newline = text.find('\n', pos);
    const size_t end = newline == std::string::npos ? text.size() : newline;
    const size_t next = newline == std::string::npos ? text.size() : newline + 1;
    const std::string raw = text.substr(pos, end - pos);
    const std::string line = FormTrim(raw);
    if (!line.empty() && line[0] == '[') {
      section = FormTrim(line.substr(1, line.size() - 2));
      if (section == target_section) found_section = true;
    } else if (!line.empty() && line[0] != '#' && line[0] != ';') {
      const size_t equals = raw.find('=');
      if (section == target_section && equals != std::string::npos &&
          FormTrim(raw.substr(0, equals)) == target_name) {
        size_t first = raw.find_first_not_of(" \t\r", equals + 1);
        if (first == std::string::npos) {
          first = raw.size();
          if (first && raw[first - 1] == '\r') --first;
          value_begin = value_end = pos + first;
        } else {
          value_begin = pos + first;
          value_end = pos + raw.find_last_not_of(" \t\r") + 1;
        }
      }
    }
    if (section == target_section) insertion = next;
    pos = next;
  }
  std::string result = text;
  if (value_begin != std::string::npos) {
    result.replace(value_begin, value_end - value_begin, value);
    return result;
  }
  const std::string eol = text.find("\r\n") == std::string::npos ? "\n" : "\r\n";
  if (insertion == std::string::npos) insertion = result.size();
  std::string addition;
  if (insertion && result[insertion - 1] != '\n') addition += eol;
  if (!found_section) addition += "[" + target_section + "]" + eol;
  addition += target_name + " = " + value + eol;
  result.insert(insertion, addition);
  return result;
}
BALTAM_PLUGIN_FCN(HymosUIConfig) {
  try {
    Need((nrhs == 2 || nrhs == 3) && nlhs == 1,
         "hymos_ui_config expects text, Section/Key, optional new value, and one output");
    const std::string text = Text(prhs[0], "configuration must be text");
    const std::string key = Text(prhs[1], "form key must be text");
    ConfigSnapshot snapshot;
    const ConfigValidation check = ValidateAndSnapshotConfigText("couette", text, "form://couette", &snapshot);
    if (!check.valid) throw std::runtime_error(check.Message());
    const auto found = snapshot.entries.find(key);
    const std::string old = found == snapshot.entries.end() ? "" : found->second;
    std::string output = old;
    if (nrhs == 3) {
      const std::string input = Text(prhs[2], "form value must be text");
      Need(input.find_first_of("\r\n") == std::string::npos &&
           input.find('\0') == std::string::npos, "form value must occupy one line");
      const std::string value = FormTrim(input);
      output = value == old ? text : PatchFormText(text, key, value);
    }
    plhs[0] = bxCreateStringScalar(output.c_str());
  } catch (const std::exception& e) { Fail(e); }
}

BALTAM_PLUGIN_FCN(HymosSetup) {
  try {
    Need(nrhs == 1 && nlhs == 0,
         "hymos_setup expects one case input and no output");
    std::string requested_case = Text(prhs[0], "case must be text");
    // Public GUI alias only: 1D uses the established Couette configuration
    // chain internally, so no solver case ID or task ABI is introduced here.
    if (requested_case == "1D" || requested_case == "1d")
      requested_case = "couette";
    const std::string plugin_dir = PluginDirectory();
    const std::string template_path =
        plugin_dir + "/assets/couette_default_input.txt";

    ConfigSnapshot snapshot;
    ConfigValidation validation =
        ValidateAndSnapshotConfig(requested_case, template_path, &snapshot);
    if (!validation.valid) throw std::runtime_error(validation.Message());
    Need(snapshot.case_id == "couette",
         "hymos_setup currently supports only 1D parallel-plate configuration");

    CouetteConfig config;
    validation = BuildCouetteConfig(snapshot, &config);
    if (!validation.valid) throw std::runtime_error(validation.Message());

    bxArray* ui_path = bxCreateStringScalar((plugin_dir + "/ui").c_str());
    const bxArray* addpath_inputs[1] = {ui_path};
    CallHostNoOutput("addpath", 1, addpath_inputs);
    bxDestroyArray(ui_path);

    bxArray* case_arg = bxCreateStringScalar(snapshot.case_id.c_str());
    bxArray* template_arg = bxCreateStringScalar(template_path.c_str());
    const bxArray* setup_inputs[2] = {case_arg, template_arg};
    CallHostNoOutput("hymos_setup_ui", 2, setup_inputs);
    bxDestroyArray(template_arg);
    bxDestroyArray(case_arg);
  } catch (const std::exception& e) {
    Fail(e);
  }
}

BALTAM_PLUGIN_FCN(HymosInfo) {
  try {
    Need(nrhs==0 && nlhs<=1,"hymos_info expects no inputs and at most one output");
    const char* names[]={"version","supported_cases","verified_plugin_paths","core_capabilities"};
    bxArray* out=bxCreateStructMatrix(1,1,4,names);
    bxSetField(out,0,"version",bxCreateStringScalar("HyMoS 1.0.0"));
    bxSetField(out,0,"supported_cases",bxCreateStringScalar("couette"));
    bxSetField(out,0,"verified_plugin_paths",bxCreateStringScalar("1D parallel-plate BGK: Couette, Fourier and Couette-Fourier through configurable wall temperatures and tangential velocities; ORDER=3..49, single/NMG, FIM=1/2/3, uniform mesh, n_thread=1..256"));
    bxSetField(out,0,"core_capabilities",bxCreateStringScalar("couette, shockstructure, cavity; FIM-1/2/3; single/NMG"));
    if(nlhs) plhs[0]=out; else bxDestroyArray(out);
  } catch(const std::exception& e) { Fail(e); }
}
BALTAM_PLUGIN_FCN(HymosValidate) {
  try {
    Need(nrhs==2 && nlhs==1,"hymos_validate expects case, configfile and one output");
    ConfigSnapshot snapshot; ConfigValidation v=ValidateAndSnapshotConfig(Text(prhs[0],"case must be text"),Text(prhs[1],"configfile must be text"),&snapshot);
    if (v.valid && snapshot.case_id == "couette") {
      CouetteConfig config; v=BuildCouetteConfig(snapshot,&config);
    } else if (v.valid) {
      v.valid=false; v.errors.push_back("case is not yet exposed by the plugin");
    }
    const char* names[]={"valid","case","message"}; bxArray* out=bxCreateStructMatrix(1,1,3,names);
    bxSetField(out,0,"valid",bxCreateLogicalScalar(v.valid));
    bxSetField(out,0,"case",bxCreateStringScalar(v.canonical_case_id.c_str()));
    bxSetField(out,0,"message",bxCreateStringScalar(v.Message().c_str())); plhs[0]=out;
  } catch(const std::exception& e) { Fail(e); }
}
BALTAM_PLUGIN_FCN(HymosValidateText) {
  try {
    Need(nrhs == 2 && nlhs == 1,
         "hymos_validate_text expects case, configuration text and one output");
    ConfigSnapshot snapshot;
    ConfigValidation v = ValidateAndSnapshotConfigText(
        Text(prhs[0], "case must be text"),
        Text(prhs[1], "configuration must be text"),
        "editor://couette", &snapshot);
    if (v.valid && snapshot.case_id == "couette") {
      CouetteConfig config;
      v = BuildCouetteConfig(snapshot, &config);
    } else if (v.valid) {
      v.valid = false;
      v.errors.push_back("case is not yet exposed by the plugin");
    }
    const char* names[] = {"valid", "case", "message"};
    bxArray* out = bxCreateStructMatrix(1, 1, 3, names);
    bxSetField(out, 0, "valid", bxCreateLogicalScalar(v.valid));
    bxSetField(out, 0, "case",
               bxCreateStringScalar(v.canonical_case_id.c_str()));
    bxSetField(out, 0, "message", bxCreateStringScalar(v.Message().c_str()));
    plhs[0] = out;
  } catch (const std::exception& e) {
    Fail(e);
  }
}
BALTAM_PLUGIN_FCN(HymosSubmit) {
  try {
    Need(nrhs==2 && nlhs==1,"hymos_submit expects case, configfile and one task output");
    ConfigSnapshot snapshot;
    ConfigValidation validation=ValidateAndSnapshotConfig(
        Text(prhs[0],"case must be text"),
        Text(prhs[1],"configfile must be text"),&snapshot);
    if(!validation.valid) throw std::runtime_error(validation.Message());
    if(snapshot.case_id!="couette")
      throw std::runtime_error("only verified couette is available");
    plhs[0]=CreateTaskObject(snapshot);
  } catch(const std::exception& e) { Fail(e); }
}
BALTAM_PLUGIN_FCN(HymosSubmitText) {
  try {
    Need(nrhs==2 && nlhs==1,
         "hymos_submit_text expects case, configuration text and one task output");
    ConfigSnapshot snapshot;
    ConfigValidation validation=ValidateAndSnapshotConfigText(
        Text(prhs[0],"case must be text"),
        Text(prhs[1],"configuration must be text"),
        "editor://couette",&snapshot);
    if(!validation.valid) throw std::runtime_error(validation.Message());
    if(snapshot.case_id!="couette")
      throw std::runtime_error("only verified couette is available");
    plhs[0]=CreateTaskObject(snapshot);
  } catch(const std::exception& e) { Fail(e); }
}
BALTAM_PLUGIN_FCN(HymosStatus) {
  try { Need(nrhs==1 && nlhs==1,"hymos_status expects task and one output"); plhs[0]=Status(GetTask(prhs[0])->value); }
  catch(const std::exception& e) { Fail(e); }
}
BALTAM_PLUGIN_FCN(HymosWait) {
  try {
    Need(nrhs==2 && nlhs==1,"hymos_wait expects task, timeout and one output");
    Need(bxIsRealDouble(prhs[1])&&bxIsScalar(prhs[1]),"timeout must be a real scalar");
    double seconds=bxGetScalar(prhs[1]); Need(std::isfinite(seconds)&&seconds>=0&&seconds<=3600,"timeout must be in [0,3600] seconds");
    Task* task=GetTask(prhs[0]); task->value->WaitFor(std::chrono::milliseconds((long long)(seconds*1000)));
    plhs[0]=Status(task->value);
  } catch(const std::exception& e) { Fail(e); }
}
BALTAM_PLUGIN_FCN(HymosStop) {
  try {
    Need(nrhs==1 && nlhs<=1,"hymos_stop expects task and at most one output");
    std::shared_ptr<CouetteTask> task=GetTask(prhs[0])->value; task->RequestStop();
    if(nlhs) plhs[0]=Status(task);
  } catch(const std::exception& e) { Fail(e); }
}
BALTAM_PLUGIN_FCN(HymosResult) {
  try {
    Need(nrhs==1 && nlhs==1,"hymos_result expects task and one output");
    std::shared_ptr<CouetteTask> task=GetTask(prhs[0])->value; CouetteStatus state=task->status();
    Need(Done(state),"task is still running; use hymos_status or hymos_wait");
    CouetteResult result=task->GetResult();
    if(state==CouetteStatus::InvalidConfig||state==CouetteStatus::Busy||state==CouetteStatus::Failed) throw std::runtime_error(result.message);
    plhs[0]=Result(result, task->progress());
  } catch(const std::exception& e) { Fail(e); }
}
BALTAM_PLUGIN_FCN(HymosExport) {
  try {
    Need(nrhs == 1 && nlhs <= 1,
         "hymos_export expects one task and at most one output");
    bxArray* summary = ExportResultFiles(GetTask(prhs[0])->value);
    if (nlhs) plhs[0] = summary; else bxDestroyArray(summary);
  } catch (const std::exception& e) {
    Fail(e);
  }
}
} // namespace

extern "C" BEX_EXPORT int bxPluginInitLib(void*) { return 0; }
extern "C" BEX_EXPORT int bxPluginInit(int,const bxArray*[]) {
  if (g_task_struct_registered) return 0;
  g_task_struct_id = bxRegisterCStruct("HyMoS", CopyTask, DeleteTask);
  bxSetCStructName(g_task_struct_id, "HyMoS task");
  g_task_struct_registered = true;
  return 0;
}
extern "C" BEX_EXPORT int bxPluginFini(void) {
  std::vector<std::shared_ptr<CouetteTask> > tasks;
  {
    std::lock_guard<std::mutex> lock(g_tasks_mutex);
    for (std::vector<std::weak_ptr<CouetteTask> >::const_iterator it = g_tasks.begin();
         it != g_tasks.end(); ++it) {
      std::shared_ptr<CouetteTask> task = it->lock();
      if (task) tasks.push_back(task);
    }
    g_tasks.clear();
  }
  for (size_t i = 0; i < tasks.size(); ++i) tasks[i]->RequestStop();
  for (size_t i = 0; i < tasks.size(); ++i) tasks[i]->Wait();
  return 0;
}
extern "C" BEX_EXPORT bexfun_info_t* bxPluginFunctions(void) {
  static bexfun_info_t f[]={
    {"hymos_info",HymosInfo,"HyMoS plugin capabilities"},
    {"hymos_run",HymosRun,"Submit applied Couette parameters and monitor in the foreground"},
    {"hymos_ui_config",HymosUIConfig,"Internal lossless form/text bridge"},
    {"hymos_ui_file",HymosUIFile,"Internal editor file adapter"},
    {"hymos_setup",HymosSetup,"Open the Couette INI configuration editor"},
    {"hymos_validate",HymosValidate,"Validate legacy HyMoS TXT configuration"},
    {"hymos_validate_text",HymosValidateText,"Internal GUI validation bridge"},
    {"hymos_submit",HymosSubmit,"Submit a non-blocking HyMoS task"},
    {"hymos_submit_text",HymosSubmitText,"Internal GUI task submission bridge"},
    {"hymos_status",HymosStatus,"Read task status"},
    {"hymos_wait",HymosWait,"Bounded task wait"},
    {"hymos_stop",HymosStop,"Request a safe task stop"},
    {"hymos_result",HymosResult,"Read structured in-memory result"},
    {"hymos_export",HymosExport,"Export Couette result files to the task workspace"},
    {NULL,NULL,NULL}};
  return f;
}
