#pragma once

#include <gkyl_basis.h>

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define APP_ARGS_DEFAULT_FILE_NAME "name-not-set"

#define APP_ARGS_CHOOSE(val, def) ((val) > 0 ? (val) : (def))

struct gkyl_app_args {
  bool use_gpu; // should this be run on GPU?
  bool step_mode; // run for fixed number of steps? (for valgrind/cuda-memcheck)
  bool trace_mem; // should we trace memory allocation/deallocations?
  int num_steps; // number of steps
  int num_threads; // number of threads
  int xcells[3]; // configuration space cells
  int vcells[3]; // velocity space cells
  char file_name[1024]; // name of input file
  enum gkyl_basis_type basis_type; // type of basis functions to use
};

static int
get_basis_type(const char *nm)
{
  if (strcmp(nm, "ms") == 0) {
    return GKYL_BASIS_MODAL_SERENDIPITY;
  }
  else if (strcmp(nm, "mt") == 0) {
    return GKYL_BASIS_MODAL_TENSOR;
  }
  return -1;
}

static struct gkyl_app_args
parse_app_args(int argc, char **argv)
{
  bool use_gpu = false;
  bool step_mode = false;
  bool trace_mem = false;
  int num_steps = INT_MAX;
  int num_threads = 1; // by default use only 1 thread

  struct gkyl_app_args args = {
    .xcells = { 0 },
    .vcells = { 0 }
  };

  strcpy(args.file_name, APP_ARGS_DEFAULT_FILE_NAME); // default
  args.basis_type = GKYL_BASIS_MODAL_SERENDIPITY;

  int c;
  while ((c = getopt(argc, argv, "+hgmt:s:i:b:x:y:z:u:v:w:")) != -1) {
    switch (c)
    {
      case 'h':
        printf("Usage: <app_name> -g -m -s nsteps -t nthreads -i inp -b [ms|mt] -x NX -y NY -z NZ -u VX -v VY -w VZ\n");
        printf(" All flags and parameters are optional.\n");
        printf(" -g     Run on GPUs if GPUs are present and code built for GPUs\n");
        printf(" -sN    Only run N steps of simulation\n");
        printf(" -tN    Use N threads (when available)\n");
        printf(" -b     Basis function to use (ms: Modal serendipity; mt: Modal tensor-product)\n");
        printf("        (Ignored for finite-volume solvers)\n");
        printf(" -m     Turn on memory allocation/deallocation tracing\n");
        printf("\n");
        printf(" Grid resolution in configuration space:\n");
        printf(" -xNX -yNY -zNZ\n");
        printf(" Grid resolution in velocity space:\n");
        printf(" -uVX -vVY -wVZ\n");
        exit(-1);
        break;

      case 'g':
        use_gpu = true;
        break;

      case 'm':
        trace_mem = true;
        break;
      
      case 's':
        step_mode = true;
        num_steps = atoi(optarg);
        break;

      case 't':
        num_threads = atoi(optarg);
        break;

      case 'x':
        args.xcells[0] = atoi(optarg);
        break;
        
      case 'y':
        args.xcells[1] = atoi(optarg);
        break;
        
      case 'z':
        args.xcells[2] = atoi(optarg);
        break;

      case 'u':
        args.vcells[0] = atoi(optarg);
        break;
        
      case 'v':
        args.vcells[1] = atoi(optarg);
        break;
        
      case 'w':
        args.vcells[2] = atoi(optarg);
        break;        

      case 'i':
        strcpy(args.file_name, optarg);
        break;

      case 'b':
        args.basis_type = get_basis_type(optarg);
        assert(args.basis_type != -1);
        break;

      case '?':
        break;
    }
  }
  
  args.use_gpu = use_gpu;
  args.trace_mem = trace_mem;
  args.step_mode = step_mode;
  args.num_steps = num_steps;
  args.num_threads = num_threads;

  return args;
}
