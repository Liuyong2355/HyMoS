#ifndef HYMOS_COMMON_CONTROLLER_H
#define HYMOS_COMMON_CONTROLLER_H

#ifndef PARALLEL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct Controller {
  int n_entry;
  struct ControllerEntry {
    char name[16];
    void (*function)(const char*);
    char description[64];
  } entry[256];
};

void resetController();
void registerController(const char* name, void (*function)(const char*), const char* description);
void clearControlRequests();
void controllerScript(const char* filename);
void getControl();

#endif

#endif
