#pragma once

#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

struct gkyl_app_args {
  bool use_gpu; // should this be run on GPU?
  bool step_mode; // run for fixed number of steps? (for valgrind/cuda-memcheck)
  int num_steps; // number of steps
  char file_name[1024]; // name of input file
};

static struct gkyl_app_args
get_parse_app_args(int argc, char **argv)
{
  bool use_gpu = false;
  bool step_mode = false;
  int num_steps = INT_MAX;

  struct gkyl_app_args args;

  int c;
  while ((c = getopt(argc, argv, "+hgs:i:")) != -1) {
    switch (c)
    {
      case 'h':
        printf("Usage: <app_name> -g -s num_steps -i inp_file \n");
        exit(-1);
        break;

      case 'g':
        use_gpu = true;
        break;        
      
      case 's':
        step_mode = true;
        num_steps = atoi(optarg);
        break;

      case 'i':
        strcpy(args.file_name, optarg);
        break;

      case '?':
        break;
    }
  }
  
  args.use_gpu = use_gpu;
  args.step_mode = step_mode;
  args.num_steps = num_steps;

  return args;
}
