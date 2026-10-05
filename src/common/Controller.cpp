#ifndef PARALLEL

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "Controller.h"

struct Controller controller;

void resetController() {
  memset(&controller, 0, sizeof(controller));
}

void registerController(const char* name, void (*function)(const char*), const char* description) {
  snprintf(controller.entry[controller.n_entry].name, sizeof(controller.entry[controller.n_entry].name), "%s", name);
  controller.entry[controller.n_entry].function = function;
  snprintf(controller.entry[controller.n_entry].description,
           sizeof(controller.entry[controller.n_entry].description),
           "%s",
           description);
  controller.n_entry++;
}

void clearControlRequests() {
  DIR* dir = opendir(".");
  if (dir == NULL) {
    return;
  }

  struct dirent* entry = NULL;
  while ((entry = readdir(dir)) != NULL) {
    if (strncmp(entry->d_name, ".control.", 9) == 0) {
      unlink(entry->d_name);
    }
  }

  closedir(dir);
}

int getOneControl(const int index) {
  int control = 0;
  int dummy = 0;
  size_t n = 0;
  char* arguments = NULL;
  char file[128];
  FILE* fp = NULL;
  sprintf(file, ".control.%d", index);
  if ((fp = fopen(file, "r")) == NULL) {
    return 0;
  }

  dummy = fscanf(fp, "%d", &control);
  dummy = getline(&arguments, &n, fp);
  fclose(fp);
  (void)dummy;
  if (control >= 0 && control < controller.n_entry) {
    (*controller.entry[control].function)(arguments);
  }
  if (arguments) {
    free(arguments);
  }

  unlink(file);
  return 1;
}

void getControl() {
  int index = 1;
  do {
    if (!getOneControl(index++)) {
      break;
    }
  } while (1);
}

void controllerScript(const char* file) {
  int i = 0;
  FILE* fp = fopen(file, "w");
  if (fp == NULL) {
    return;
  }

  fprintf(fp, "#! /bin/bash\n");
  fprintf(fp, "\n");
  fprintf(fp, "RETVAL=0\n");
  fprintf(fp, "for ((i = 1;i > 0;i = i+1)); do \n");
  fprintf(fp, "  if ! [ -f .control.$i ]; then\n");
  fprintf(fp, "    file=\".control.$i\" ;\n");
  fprintf(fp, "    echo \"This is the No.\" $i \" unmanaged requests.\" ;\n");
  fprintf(fp, "    i=-1 ;\n");
  fprintf(fp, "  fi\n");
  fprintf(fp, "done\n\n");
  fprintf(fp, "case \"$1\" in\n");

  for (i = 0; i < controller.n_entry; i++) {
    fprintf(fp, "  %s)\n", controller.entry[i].name);
    fprintf(fp, "    echo \"%d $2 $3 $4 $5\" >$file\n", i);
    fprintf(fp, "    ;;\n");
  }

  fprintf(fp, "  *)\n");
  fprintf(fp, "    echo $\"Unmatched argument, valid arguments are:\"\n");
  for (i = 0; i < controller.n_entry; i++) {
    fprintf(fp, "    echo \"\t%s: %s\"\n", controller.entry[i].name, controller.entry[i].description);
  }
  fprintf(fp, "    RETVAL=1\n");
  fprintf(fp, "esac\n");
  fprintf(fp, "\n");
  fprintf(fp, "exit $RETVAL\n");
  fclose(fp);

  chmod(file, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
  const char* display_file = getenv("HYMOS_CONTROL_DISPLAY");
  if (display_file == NULL || display_file[0] == '\0') {
    display_file = file;
  }
  printf("\nUsing %s to manage the job.\n", display_file);
}

#endif
