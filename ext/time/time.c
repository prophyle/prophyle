/* `time' utility to display resource usage of processes.
   Copyright (C) 1990, 91, 92, 93, 96 Free Software Foundation, Inc.

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

/* Originally written by David Keppel <pardo@cs.washington.edu>.
   Heavily modified by David MacKenzie <djm@gnu.ai.mit.edu>.  */

extern const char *version_string;

#include "wait.h"
#include <stdio.h>
#include <sys/param.h>		/* For getpagesize, maybe.  */
#include <signal.h>
#include <errno.h>
#ifndef errno
extern int errno;
#endif
#include <getopt.h>
#include "port.h"
#include "resuse.h"
#include "getpagesize.h"

void error PARAMS((int status, int errnum, char *message, ...));

static void usage PARAMS((FILE *, int));

/* A Pointer to a signal handler.  */
typedef RETSIGTYPE (*sighandler) ();

/* msec = milliseconds = 1/1,000 (1*10e-3) second.
   usec = microseconds = 1/1,000,000 (1*10e-6) second.  */

/* Systems known to fill in the average resident set size fields:
   SunOS 4.1.3 (m68k and sparc)
   Mt. Xinu 4.3BSD on HP9000/300 (m68k)
   Ultrix 4.4 (mips)
   IBM ACIS 4.3BSD (rt)
   Sony NEWS-OS 4.1C (m68k)

   Systems known to not fill them in:
   OSF/1 1.3 (alpha)
   BSD/386 1.1 (anything derived from NET-2)
   NetBSD 1.0 (4.4BSD-derived)
   Irix 5.2 (R4000)
   Solaris 2.3
   Linux 1.0

   It doesn't matter how many clock ticks/second there are on
   systems that don't fill in those fields.

   If the avgresident (%t) we print is less than a power of 2 away from
   the maxresident (%M), then we likely are using the right number.
   Another good check is comparing the average text size with the
   output of `size' on the executable.

   According to the SunOS manual, there are 50 ticks/sec on the sun3
   and 100 on the sun4.

   Some manuals have an apparent error, claiming that units for average
   sizes are kb*sec.  Judging by the contents of `struct rusage', it
   looks like it should be kb*ticks, like on SunOS.  Ticks/sec seems
   to be (empirically):
   50 Mt. Xinu
   250 Ultrix (mips)
   50 ACIS
   100 NEWS-OS

   sysconf(_SC_CLK_TCK) is *unrelated*.  */

#if defined(sun3) || defined(hp300) || defined(ibm032)
#define TICKS_PER_SEC 50
#endif
#if defined(mips)
#define TICKS_PER_SEC 250
#endif
#ifndef TICKS_PER_SEC
#define TICKS_PER_SEC 100
#endif

/* The number of milliseconds in one `tick' used by the `rusage' structure.  */
#define MSEC_PER_TICK (1000 / TICKS_PER_SEC)

/* Return the number of clock ticks that occur in M milliseconds.  */
#define MSEC_TO_TICKS(m) ((m) / MSEC_PER_TICK)

typedef enum {false, true} boolean;

#define UL unsigned long

/* The default output format.  */
static const char *const default_format =
"%Uuser %Ssystem %Eelapsed %PCPU (%Xavgtext+%Davgdata %Mmaxresident)k\n\
%Iinputs+%Ooutputs (%Fmajor+%Rminor)pagefaults %Wswaps";

/* The output format for the -p option .*/
static const char *const posix_format = "real %e\nuser %U\nsys %S";

/* Format string for printing all statistics verbosely.
   Keep this output to 24 lines so users on terminals can see it all.

   The format string is used two ways: as a format string, and in
   verbose mode, to document all the possible formatting possiblities.
   When `longstats' is used as a format string, it has to be put into
   one contiguous string (e.g., into a `char[]').  We could alternatively
   store it as a `char *' and convert it into a `*char[]' when we need
   it as documentation, but most compilers choke on very long strings.  */

static const char *const longstats[] =
{
  "\tCommand being timed: \"%C\"\n",
  "\tUser time (seconds): %U\n",
  "\tSystem time (seconds): %S\n",
  "\tPercent of CPU this job got: %P\n",
  "\tElapsed (wall clock) time (h:mm:ss or m:ss): %E\n",
  "\tAverage shared text size (kbytes): %X\n",
  "\tAverage unshared data size (kbytes): %D\n",
  "\tAverage stack size (kbytes): %p\n",
  "\tAverage total size (kbytes): %K\n",
  "\tMaximum resident set size (kbytes): %M\n",
  "\tAverage resident set size (kbytes): %t\n",
  "\tMajor (requiring I/O) page faults: %F\n",
  "\tMinor (reclaiming a frame) page faults: %R\n",
  "\tVoluntary context switches: %w\n",
  "\tInvoluntary context switches: %c\n",
  "\tSwaps: %W\n",
  "\tFile system inputs: %I\n",
  "\tFile system outputs: %O\n",
  "\tSocket messages sent: %s\n",
  "\tSocket messages received: %r\n",
  "\tSignals delivered: %k\n",
  "\tPage size (bytes): %Z\n",
  "\tExit status: %x",
  NULL
};

/* If true, show an English description next to each statistic.  */
static boolean verbose;

/* Name of output file.  Only used if -o option is given.  */
static const char *outfile;

/* Output stream, stderr by default.  */
static FILE *outfp;

/* If true, append to `outfile' rather than truncating it.  */
static boolean append;

/* The output format string.  */
static const char *output_format;

/* The name this program was invoked by.  */
char *program_name;

static struct option longopts[] =
{
  {"append", no_argument, NULL, 'a'},
  {"format", required_argument, NULL, 'f'},
  {"help", no_argument, NULL, 'h'},
  {"output-file", required_argument, NULL, 'o'},
  {"portability", no_argument, NULL, 'p'},
  {"verbose", no_argument, NULL, 'v'},
  {"version", no_argument, NULL, 'V'},
  {NULL, no_argument, NULL, 0}
};

/* Print ARGV to FP, with each entry in ARGV separated by FILLER.  */

static void
fprintargv (fp, argv, filler)
     FILE *fp;
     const char *const *argv;
     const char *filler;
{
  const char *const *av;

  av = argv;
  fputs (*av, fp);
  while (*++av)
    {
      fputs (filler, fp);
      fputs (*av, fp);
    }
  if (ferror (fp))
    error (1, errno, "write error");
}

/* Return a null-terminated string containing the concatenation,
   in order, of all of the elements of ARGV.
   The '\0' at the end of each ARGV-element is not copied.
   Example:	char *argv[] = {"12", "ab", ".,"};
 		linear_argv(argv) == "12ab.,"
   Print a message and return NULL if memory allocation failed.  */

static char *
linear_argv (argv)
     const char *const *argv;
{
  const char *const *s;		/* Each string in ARGV.  */
  char *new;			/* Allocated space.  */
  char *dp;			/* Copy in to destination.  */
  const char *sp;		/* Copy from source.  */
  int size;

  /* Find length of ARGV and allocate.  */
  size = 1;
  for (s = argv; *s; ++s)
    size += strlen (*s);
  new = (char *) malloc (size);
  if (new == NULL)
    {
      fprintf (stderr, "%s: virtual memory exhausted\n", program_name);
      return NULL;
    }

  /* Copy each string in ARGV to the new string.  At the end of
     each string copy, back up `dp' so that on the next string,
     the `\0' will be overwritten.  */
  for (s = argv, sp = *s, dp = new; *s; ++s)
    {
      sp = *s;
      while ((*dp++ = *sp++) != '\0')
	/* Do nothing.  */ ;
      --dp;
    }

  return new;
}

/* Return the number of kilobytes corresponding to a number of pages PAGES.
   (Actually, we use it to convert pages*ticks into kilobytes*ticks.)

   Try to do arithmetic so that the risk of overflow errors is minimized.
   This is funky since the pagesize could be less than 1K.
   Note: Some machines express getrusage statistics in terms of K,
   others in terms of pages.  */

static unsigned long
ptok (pages)
     unsigned long pages;
{
  static unsigned long ps = 0;
  unsigned long tmp;
  static long size = LONG_MAX;

  /* Initialization.  */
  if (ps == 0)
    ps = (long) getpagesize ();

  /* Conversion.  */
  if (pages > (LONG_MAX / ps))
    {				/* Could overflow.  */
      tmp = pages / 1024;	/* Smaller first, */
      size = tmp * ps;		/* then larger.  */
    }
  else
    {				/* Could underflow.  */
      tmp = pages * ps;		/* Larger first, */
      size = tmp / 1024;	/* then smaller.  */
    }
  return size;
}

/* summarize: Report on the system use of a command.

   Copy the FMT argument to FP except that `%' sequences
   have special meaning, and `\n' and `\t' are translated into
   newline and tab, respectively, and `\\' is translated into `\'.

   The character following a `%' can be:
   (* means the tcsh time builtin also recognizes it)
   % == a literal `%'
   C == command name and arguments
*  D == average unshared data size in K (ru_idrss+ru_isrss)
*  E == elapsed real (wall clock) time in [hour:]min:sec
*  F == major page faults (required physical I/O) (ru_majflt)
*  I == file system inputs (ru_inblock)
*  K == average total mem usage (ru_idrss+ru_isrss+ru_ixrss)
*  M == maximum resident set size in K (ru_maxrss)
*  O == file system outputs (ru_oublock)
*  P == percent of CPU this job got (total cpu time / elapsed time)
*  R == minor page faults (reclaims; no physical I/O involved) (ru_minflt)
*  S == system (kernel) time (seconds) (ru_stime)
*  U == user time (seconds) (ru_utime)
*  W == times swapped out (ru_nswap)
*  X == average amount of shared text in K (ru_ixrss)
   Z == page size
*  c == involuntary context switches (ru_nivcsw)
   e == elapsed real time in seconds
*  k == signals delivered (ru_nsignals)
   p == average unshared stack size in K (ru_isrss)
*  r == socket messages received (ru_msgrcv)
*  s == socket messages sent (ru_msgsnd)
   t == average resident set size in K (ru_idrss)
*  w == voluntary context switches (ru_nvcsw)
   x == exit status of command

   Various memory usages are found by converting from page-seconds
   to kbytes by multiplying by the page size, dividing by 1024,
   and dividing by elapsed real time.

   FP is the stream to print to.
   FMT is the format string, interpreted as described above.
   COMMAND is the command and args that are being summarized.
   RESP is resource information on the command.  */

static void
summarize (fp, fmt, command, resp)
     FILE *fp;
     const char *fmt;
     const char **command;
     RESUSE *resp;
{
  unsigned long r;		/* Elapsed real milliseconds.  */
  unsigned long v;		/* Elapsed virtual (CPU) milliseconds.  */

  if (WIFSTOPPED (resp->waitstatus))
    fprintf (fp, "Command stopped by signal %d\n",
	     WSTOPSIG (resp->waitstatus));
  else if (WIFSIGNALED (resp->waitstatus))
    fprintf (fp, "Command terminated by signal %d\n",
	     WTERMSIG (resp->waitstatus));
  else if (WIFEXITED (resp->waitstatus) && WEXITSTATUS (resp->waitstatus))
    fprintf (fp, "Command exited with non-zero status %d\n",
	     WEXITSTATUS (resp->waitstatus));

  /* Convert all times to milliseconds.  Occasionally, one of these values
     comes out as zero.  Dividing by zero causes problems, so we first
     check the time value.  If it is zero, then we take `evasive action'
     instead of calculating a value.  */

  r = resp->elapsed.tv_sec * 1000 + resp->elapsed.tv_usec / 1000;

  v = resp->ru.ru_utime.tv_sec * 1000 + resp->ru.ru_utime.TV_MSEC +
    resp->ru.ru_stime.tv_sec * 1000 + resp->ru.ru_stime.TV_MSEC;

  while (*fmt)
    {
      switch (*fmt)
	{
	case '%':
	  switch (*++fmt)
	    {
	    case '%':		/* Literal '%'.  */
	      putc ('%', fp);
	      break;
	    case 'C':		/* The command that got timed.  */
	      fprintargv (fp, command, " ");
	      break;
	    case 'D':		/* Average unshared data size.  */
	      fprintf (fp, "%lu",
		       MSEC_TO_TICKS (v) == 0 ? 0 :
		       ptok ((UL) resp->ru.ru_idrss) / MSEC_TO_TICKS (v) +
		       ptok ((UL) resp->ru.ru_isrss) / MSEC_TO_TICKS (v));
	      break;
	    case 'E':		/* Elapsed real (wall clock) time.  */
	      if (resp->elapsed.tv_sec >= 3600)	/* One hour -> h:m:s.  */
		fprintf (fp, "%ld:%02ld:%02ld",
			 resp->elapsed.tv_sec / 3600,
			 (resp->elapsed.tv_sec % 3600) / 60,
			 resp->elapsed.tv_sec % 60);
	      else
		fprintf (fp, "%ld:%02ld.%02ld",	/* -> m:s.  */
			 resp->elapsed.tv_sec / 60,
			 resp->elapsed.tv_sec % 60,
			 resp->elapsed.tv_usec / 10000);
	      break;
	    case 'F':		/* Major page faults.  */
	      fprintf (fp, "%ld", resp->ru.ru_majflt);
	      break;
	    case 'I':		/* Inputs.  */
	      fprintf (fp, "%ld", resp->ru.ru_inblock);
	      break;
	    case 'K':		/* Average mem usage == data+stack+text.  */
	      fprintf (fp, "%lu",
		       MSEC_TO_TICKS (v) == 0 ? 0 :
		       ptok ((UL) resp->ru.ru_idrss) / MSEC_TO_TICKS (v) +
		       ptok ((UL) resp->ru.ru_isrss) / MSEC_TO_TICKS (v) +
		       ptok ((UL) resp->ru.ru_ixrss) / MSEC_TO_TICKS (v));
	      break;
	    case 'M':		/* Maximum resident set size.  */
	      fprintf (fp, "%lu", ptok ((UL) resp->ru.ru_maxrss));
	      break;
	    case 'O':		/* Outputs.  */
	      fprintf (fp, "%ld", resp->ru.ru_oublock);
	      break;
	    case 'P':		/* Percent of CPU this job got.  */
	      /* % cpu is (total cpu time)/(elapsed time).  */
	      if (r > 0)
		fprintf (fp, "%lu%%", (v * 100 / r));
	      else
		fprintf (fp, "?%%");
	      break;
	    case 'R':		/* Minor page faults (reclaims).  */
	      fprintf (fp, "%ld", resp->ru.ru_minflt);
	      break;
	    case 'S':		/* System time.  */
	      fprintf (fp, "%ld.%02ld",
		       resp->ru.ru_stime.tv_sec,
		       resp->ru.ru_stime.TV_MSEC / 10);
	      break;
	    case 'U':		/* User time.  */
	      fprintf (fp, "%ld.%02ld",
		       resp->ru.ru_utime.tv_sec,
		       resp->ru.ru_utime.TV_MSEC / 10);
	      break;
	    case 'W':		/* Times swapped out.  */
	      fprintf (fp, "%ld", resp->ru.ru_nswap);
	      break;
	    case 'X':		/* Average shared text size.  */
	      fprintf (fp, "%lu",
		       MSEC_TO_TICKS (v) == 0 ? 0 :
		       ptok ((UL) resp->ru.ru_ixrss) / MSEC_TO_TICKS (v));
	      break;
	    case 'Z':		/* Page size.  */
	      fprintf (fp, "%d", getpagesize ());
	      break;
	    case 'c':		/* Involuntary context switches.  */
	      fprintf (fp, "%ld", resp->ru.ru_nivcsw);
	      break;
	    case 'e':		/* Elapsed real time in seconds.  */
	      fprintf (fp, "%ld.%02ld",
		       resp->elapsed.tv_sec,
		       resp->elapsed.tv_usec / 10000);
	      break;
	    case 'k':		/* Signals delivered.  */
	      fprintf (fp, "%ld", resp->ru.ru_nsignals);
	      break;
	    case 'p':		/* Average stack segment.  */
	      fprintf (fp, "%lu",
		       MSEC_TO_TICKS (v) == 0 ? 0 :
		       ptok ((UL) resp->ru.ru_isrss) / MSEC_TO_TICKS (v));
	      break;
	    case 'r':		/* Incoming socket messages received.  */
	      fprintf (fp, "%ld", resp->ru.ru_msgrcv);
	      break;
	    case 's':		/* Outgoing socket messages sent.  */
	      fprintf (fp, "%ld", resp->ru.ru_msgsnd);
	      break;
	    case 't':		/* Average resident set size.  */
	      fprintf (fp, "%lu",
		       MSEC_TO_TICKS (v) == 0 ? 0 :
		       ptok ((UL) resp->ru.ru_idrss) / MSEC_TO_TICKS (v));
	      break;
	    case 'w':		/* Voluntary context switches.  */
	      fprintf (fp, "%ld", resp->ru.ru_nvcsw);
	      break;
	    case 'x':		/* Exit status.  */
	      fprintf (fp, "%d", WEXITSTATUS (resp->waitstatus));
	      break;
	    case '\0':
	      putc ('?', fp);
	      return;
	    default:
	      putc ('?', fp);
	      putc (*fmt, fp);
	    }
	  ++fmt;
	  break;

	case '\\':		/* Format escape.  */
	  switch (*++fmt)
	    {
	    case 't':
	      putc ('\t', fp);
	      break;
	    case 'n':
	      putc ('\n', fp);
	      break;
	    case '\\':
	      putc ('\\', fp);
	      break;
	    default:
	      putc ('?', fp);
	      putc ('\\', fp);
	      putc (*fmt, fp);
	    }
	  ++fmt;
	  break;

	default:
	  putc (*fmt++, fp);
	}

      if (ferror (fp))
	error (1, errno, "write error");
    }
  putc ('\n', fp);

  if (ferror (fp))
    error (1, errno, "write error");
}

/* Initialize the options and parse the command line arguments.
   Also note the position in ARGV where the command to time starts.

   By default, output is to stderr.

   ARGV is the array of command line arguments.
   ARGC is the number of command line arguments.

   Return the command line to run and gather statistics on.  */

static const char **
getargs (argc, argv)
     int argc;
     char **argv;
{
  int optc;
  char *format;			/* Format found in environment.  */

  /* Initialize the option flags.  */
  verbose = false;
  outfile = NULL;
  outfp = stderr;
  append = false;
  output_format = default_format;
  program_name = argv[0];

  /* Set the format string from the environment.  Do this before checking
     the args so that we won't clobber a user-specified format.  */
  format = getenv ("TIME");
  if (format)
    output_format = format;

  while ((optc = getopt_long (argc, argv, "+af:o:pvV", longopts, (int *) 0))
	 != EOF)
    {
      switch (optc)
	{
	case 'a':
	  append = true;
	  break;
	case 'f':
	  output_format = optarg;
	  break;
	case 'h':
	  usage (stdout, 0);
	case 'o':
	  outfile = optarg;
	  break;
	case 'p':
	  output_format = posix_format;
	  break;
	case 'v':
	  verbose = true;
	  break;
	case 'V':
	  fprintf (stderr, "%s\n", version_string);
	  exit (0);
	default:
	  usage (stderr, 1);
	}
    }

  if (optind == argc)
    usage (stderr, 1);

  if (outfile)
    {
      if (append)
	outfp = fopen (outfile, "a");
      else
	outfp = fopen (outfile, "w");
      if (outfp == NULL)
	error (1, errno, "%s", outfile);
    }

  /* If the user specified verbose output, we need to convert
     `longstats' to a `char *'.  */
  if (verbose)
    {
      output_format = (const char *) linear_argv (longstats);
      if (output_format == NULL)
	exit (1);		/* Out of memory.  */
    }

  return (const char **) &argv[optind];
}

/* Run command CMD and return statistics on it.
   Put the statistics in *RESP.  */

static void
run_command (cmd, resp)
     char *const *cmd;
     RESUSE *resp;
{
  pid_t pid;			/* Pid of child.  */
  sighandler interrupt_signal, quit_signal;

  resuse_start (resp);

  pid = fork ();		/* Run CMD as child process.  */
  if (pid < 0)
    error (1, errno, "cannot fork");
  else if (pid == 0)
    {				/* If child.  */
      /* Don't cast execvp arguments; that causes errors on some systems,
	 versus merely warnings if the cast is left off.  */
      execvp (cmd[0], cmd);
      error (0, errno, "cannot run %s", cmd[0]);
      _exit (errno == ENOENT ? 127 : 126);
    }

  /* Have signals kill the child but not self (if possible).  */
  interrupt_signal = signal (SIGINT, SIG_IGN);
  quit_signal = signal (SIGQUIT, SIG_IGN);

  if (resuse_end (pid, resp) == 0)
    error (1, errno, "error waiting for child process");

  /* Re-enable signals.  */
  signal (SIGINT, interrupt_signal);
  signal (SIGQUIT, quit_signal);
}

void
main (argc, argv)
     int argc;
     char **argv;
{
  const char **command_line;
  RESUSE res;

  command_line = getargs (argc, argv);
  run_command (command_line, &res);
  summarize (outfp, output_format, command_line, &res);
  fflush (outfp);

  if (WIFSTOPPED (res.waitstatus))
    exit (WSTOPSIG (res.waitstatus));
  else if (WIFSIGNALED (res.waitstatus))
    exit (WTERMSIG (res.waitstatus));
  else if (WIFEXITED (res.waitstatus))
    exit (WEXITSTATUS (res.waitstatus));
}

static void
usage (stream, status)
     FILE *stream;
     int status;
{
  fprintf (stream, "\
Usage: %s [-apvV] [-f format] [-o file] [--append] [--verbose]\n\
       [--portability] [--format=format] [--output=file] [--version]\n\
       [--help] command [arg...]\n",
	   program_name);
  exit (status);
}
