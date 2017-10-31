/*  $Id: signalqq.c,v 1.3 2012/06/27 10:10:37 mjlaine Exp $ */

/* Gfortran signal handler for ctrl-c */
/* Part of the mcmc library */

/*
 * Marko Laine 2008 <marko.laine@fmi.fi>
 * Copyrights licensed under a MIT License.
 * See the accompanying LICENSE.txt file for terms.
 */

#ifdef WIN32

/* do nothing if using (mingw) gfortran in Windows */

void signalqq_(void (*handler)()){}

#else

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>

void def_handler(int sig);

typedef void (* funptr)(int);

static funptr gf_handler = def_handler;

void sig_handler(int j, siginfo_t *info, void *dummy)
{
  /* Call fortran handler given as argument to signalqq */
  gf_handler(j);

}

/* The default signal handler */
void def_handler(int sig)
{

  /* do nothing */

}

void signalqq_(void (*handler)())
{

  struct sigaction act;

  gf_handler = (*handler); /*  save routine pointer to a static variable */

  act.sa_sigaction=(*sig_handler);
  sigemptyset(&act.sa_mask);
  act.sa_flags=SA_SIGINFO;
  sigaction(SIGHUP, &act, NULL); 
  sigaction(SIGINT, &act, NULL); 
  sigaction(SIGTERM, &act, NULL); 
  sigaction(SIGTSTP, &act, NULL); 
  sigaction(SIGUSR1, &act, NULL); 
  sigaction(SIGUSR2, &act, NULL); 

}

#endif
