// Usage:
// - Used by 1D backends to expose an interactive control loop.
// - Start the process first, then send commands through the generated script.

#ifndef HYMOS_APPS_SUPPORT_CASERUNNER1D_H
#define HYMOS_APPS_SUPPORT_CASERUNNER1D_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <unistd.h>

#include "Controller.h"

namespace hymos {
namespace apps {

template <typename APP>
class CaseRunner1D {
public:
  static int Main(APP& app, const char* script_name) {
    app_ = &app;
    resetController();

    registerController("config", &ReadConfig, "read config file. Parameters: filename");
    registerController("dump", &DumpData, "dump data. Parameters: filename");
    registerController("load", &LoadData, "load data. Parameters: filename");
    registerController("order", &ResetOrder, "reset order. Parameters: order");
    registerController("maxwellian", &Maxwellian, "set maxwellian as initial value.");
    registerController("run", &Run, "run.");
    registerController("stop", &StopRun, "stop run.");
    registerController("exit", &Stop, "exit the program.");
    registerController("pause", &Pause, "pause.");
    registerController("tolerance", &ResetTol, "reset the tolerance.");
    registerController("print_config", &PrintConfig, "print current configuration.");
    registerController("mesh", &ResetMesh, "reset mesh size.");
    registerController("mesh_levels", &ResetMeshLevels, "refine/reset mesh levels.");

    clearControlRequests();
    controllerScript(script_name);

    do {
      sleep(1);
      getControl();
    } while (!app_->is_exit);

    return 0;
  }

private:
  static void ReadConfig(const char* dummy) {
    char filename[1024];
    sscanf(dummy, " %s ", filename);
    app_->read_config(filename);
  }

  static void DumpData(const char* dummy) {
    char filename[1024];
    sscanf(dummy, " %s ", filename);
    app_->dump_data(filename);
  }

  static void LoadData(const char* dummy) {
    char filename[1024];
    sscanf(dummy, " %s ", filename);
    app_->load_data(filename);
  }

  static void ResetOrder(const char* dummy) {
    app_->reset_order(atoi(dummy));
  }

  static void ResetMesh(const char* dummy) {
    const int nx = (dummy == NULL) ? -1 : atoi(dummy);
    app_->reset_mesh(nx);
  }

  static void ResetMeshLevels(const char* dummy) {
    app_->reset_mesh_levels(atoi(dummy));
  }

  static void Maxwellian(const char* /*dummy*/) {
    app_->set_maxwellian();
  }

  static void Run(const char* /*dummy*/) {
    if (app_->is_stop_run) {
      app_->run();
    }
  }

  static void StopRun(const char* /*dummy*/) {
    app_->stop_run();
  }

  static void Stop(const char* /*dummy*/) {
    app_->is_exit = true;
  }

  static void Pause(const char* /*dummy*/) {
    std::cout << "Press Enter to continue ..." << std::flush;
    getchar();
  }

  static void PrintConfig(const char* dummy) {
    app_->print_config();
    Pause(dummy);
  }

  static void ResetTol(const char* dummy) {
    double tol = 0.0;
    sscanf(dummy, "%le", &tol);
    app_->reset_tol(tol);
  }

private:
  static APP* app_;
};

template <typename APP>
APP* CaseRunner1D<APP>::app_ = NULL;

}  // namespace apps
}  // namespace hymos

#endif
