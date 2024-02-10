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
  struct radiating_state *all_states; // indecies are z-1,charge_state
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
 * @return pointer to structure of all_radiation_states
 */
struct all_radiation_states* gkyl_read_rad_fit_params();

/* Function to return the fit information for a specfied atomic number, charge state, and ne
 * @param all_radiation_states rad data: Struct containing radiation fit data
 * @param atomic_z: Z of element for desired fit information
 * @param charge_state: charge state of element for desired fit information
 * @param a, alpha, beta, gamma, V0: fit parameters to be returned 
 * @param num_densities: maximum number of densities to return fit parameters for
 * @return 1 if fit doesn't exist
 * Note: Untested for num_densities>1
 */
int gkyl_get_fit_params(const struct all_radiation_states rad_data, int atomic_z, int charge_state, double *a, double *alpha, double *beta, double *gamma, double *V0, int num_densities);

/* Function to return the fit emissivity (Lz) and temperature closest to a given input temperature
 * @param all_radiation_states rad data: Struct containing radiation fit data
 * @param atomic_z: Z of element for desired fit information
 * @param charge_state: charge state of element for desired fit information
 * @param ne: The returned te and lz are for density closest to ne.
 * @param te: returns closest fit temperature to input te
 * @param Lz: returns Lz for closest temperature to input te 
 */
int gkyl_get_fit_lz(const struct all_radiation_states rad_data, int atomic_z, int charge_state, double ne, double* te, double* Lz);

/*  Free memory of all_radiation_states
 * @param all_radiation_states rad data: Struct containing radiation fit data
 */
void gkyl_release_fit_params(struct all_radiation_states *rad_data);
