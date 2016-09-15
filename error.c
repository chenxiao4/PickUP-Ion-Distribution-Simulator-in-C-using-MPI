#include "commonincludes.h"

static void err_doit(int, int, const char *, va_list);



void err_say(const char *fmt,...){

  va_list ap;
  
  va_start(ap,fmt);
  err_doit(1,errno,fmt,ap);
  va_end(ap);
  exit(EXIT_FAILURE);
}


static void
err_doit(int errnoflag, int error, const char *fmt, va_list ap){

  char buf[MAXLINE];
  
  vsnprintf(buf,MAXLINE-1,fmt,ap);
  if (errnoflag)
    snprintf(buf+strlen(buf),MAXLINE-strlen(buf)-1,": %s",strerror(error));
  strcat(buf,"\n");
  fflush(stdout);
  fputs(buf,stderr);
  fflush(NULL);
}


void err_exit(char *str){

  fprintf(stdout,"%s\n",str);
  exit(EXIT_FAILURE);

}


void message(char *str){

  fprintf(stdout,"%s\n",str);

}


void err_check(char *str){

  char c;

  message(str);
  
  c = getchar();

}
