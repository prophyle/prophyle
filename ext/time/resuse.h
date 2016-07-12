/* resuse.h - declarations for child process resource use library
   Copyright (C) 1993, 1996 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
   02111-1307, USA.  */

#ifndef _RESUSE_H
#define _RESUSE_H 1

#if HAVE_TIMEVAL
# include <sys/time.h>
#else
/* High resolution clock structure.  */
struct timeval
{
  long tv_sec;			/* Seconds.  */
  long tv_usec;			/* Microseconds.  */
};
#endif

#if HAVE_SYS_RUSAGE_H
/* This rusage structure measures nanoseconds instead of microseconds.  */
# define TV_MSEC tv_nsec / 1000000
# include <sys/rusage.h>
#else
# define TV_MSEC tv_usec / 1000
# if HAVE_WAIT3
#  include <sys/resource.h>
# else
/* Process resource usage structure.  */
struct rusage
{
  struct timeval ru_utime;	/* User time used.  */
  struct timeval ru_stime;	/* System time used.  */
  int ru_maxrss, ru_ixrss, ru_idrss, ru_isrss,
  ru_minflt, ru_majflt, ru_nswap, ru_inblock, 
  ru_oublock, ru_msgsnd, ru_msgrcv, ru_nsignals,
  ru_nvcsw, ru_nivcsw;
};
# endif
#endif

/* Information on the resources used by a child process.  */
typedef struct
{
  int waitstatus;
  struct rusage ru;
  struct timeval start, elapsed; /* Wallclock time of process.  */
} RESUSE;

/* Prepare to measure a child process.  */
void resuse_start PARAMS ((RESUSE *resp));

/* Wait for and fill in data on child process PID.  */
int resuse_end PARAMS ((pid_t pid, RESUSE *resp));

#endif /* _RESUSE_H */
