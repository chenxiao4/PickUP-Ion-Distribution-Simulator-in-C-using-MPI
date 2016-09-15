#ifndef MENU_H
#define MENU_H

#include "parameter.h"

static const char HFILE[] = "MENU_ACT_SWICS_";
static const char EXT[] = ".DAT";

char *input_fname(int year);
void read_input(int year, PARA_t *input);
int trim_line(char *buf);

#endif
