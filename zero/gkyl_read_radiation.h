#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#define BUFFER_LEN 100
/* Data for a single radiation fit
 */
struct rad_fit_parameters{
  double electron_density;
  double A;
  double alpha;
  double beta;
  double gamma;
  double V0;
  int te_intervals;
  double *te;
  double *Lz;
};

/* Radiation data for a single charge state, may (or may not) have fits for more than
 *    1 electron density
 */
struct radiating_state{
  bool state_exists;
  int atomic_number;
  int charge_state;
  int number_of_densities;
  double *electron_densities; // The different densities for this charge state
  struct rad_fit_parameters *rad_fits; // 1 fit per density for this charge state
};

/* All of the radiation fits read.
 *
 */
struct all_radiation_states{
  int max_atomic_number;
  struct radiating_state *all_states;
};

/* Function to read a line with two numbers in the format:
 * STRING 'DELIM' #1 'DELIM' STRING 'DELIM' #2 
 * DELIM can be any of: ,:;=
 */ 
static inline void read_two_numbers(FILE *fptr, int *num1, int *num2){
  char str[BUFFER_LEN];
  char delim[5]="=,;:";
  if(fgets(str,BUFFER_LEN,fptr)!=NULL) {
    strtok(str,delim);
    *num1=atoi(strtok(NULL,delim));
    strtok(NULL,delim);
    *num2=atoi(strtok(NULL,delim));
  }
}

/* Function to read in all the radiation fit parameters stored in "radiation_fit_params.txt"
 * 
 */
struct all_radiation_states* gkyl_read_rad_fit_params();

/* Function to return the fit information for a specfied atomic number, charge state, and ne
 *
 * Note: Untested for num_densities>1
 */
int gkyl_get_fit_params(const struct all_radiation_states rad_data, int atomic_z, int charge_state, double *a, double *alpha, double *beta, double *gamma, double *V0, int num_densities);

void gkyl_release_fit_params(struct all_radiation_states *rad_data);
