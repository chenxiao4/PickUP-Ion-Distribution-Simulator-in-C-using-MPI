#ifndef HEADER_H
#define HEADER_H

#define _POSIX_C_SOURCE 200809L
#define _XOPEN_SOURCE 700
#define _ISOC99_SOURCE
#define __EXTENSIONS__

//header file lib
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <tgmath.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>


//define minmax()
#define min(a,b) ((a) < (b)? (a) : (b))
#define max(a,b) ((a) > (b)? (a) : (b))

#define MAXLINE 4096



#endif
