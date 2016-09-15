#ifndef ERROR_H
#define ERROR_H

//#include "header.h"

void err_say(const char *,...) __attribute__((noreturn));
void err_exit(char *str);
void message(char *str);
void err_check(char *str);

#endif
