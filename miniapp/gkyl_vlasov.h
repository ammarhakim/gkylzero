#pragma once

#ifndef GKYL_MAX_SPECIES
# define GKYL_MAX_SPECIES 2
#endif

// Parameters for Vlasov species
struct gkyl_vlasov_species {
    char name[128]; // Species name
    double charge, mass; // Charge and mass
    double lower[3], upper[3]; // Lower, upper bounds of velocity-space
    int cells[3]; // Velocity-space cells

    int evolve; // Evolve species? 1-yes, 0-no

    void *ctx; // Context for initial condition init function
    // Pointer to initialization function
    void (*init)(double t, const double *xn, double *fout, void *ctx);

    int num_diag_moments; // Number of diagnostic moments
    char diag_moments[16][16]; // List of diagnostic moments
};

// Parameter for EM field
struct gkyl_em_field {
    double epsilon0, mu0;

    int evolve; // Evolve field? 1-yes, 0-no

    void *ctx; // Context for initial condition init function
    // Pointer to initialization function
    void (*init)(double t, const double *xn, double *fout, void *ctx);
};

// Top-level app parameters
struct gkyl_vm {
    char name[128]; // Name of app: used as output prefix

    int cdim, vdim; // Conf, velocity space dimensions
    double lower[3], upper[3]; // Lower, upper bounds of config-space
    int cells[3]; // Config-space cells
    int poly_order; // Polynomial order

    int num_species; // Number of species
    struct gkyl_vlasov_species species[GKYL_MAX_SPECIES]; // Species objects
    struct gkyl_em_field field; // Field object
};

// Object representing Vlasov mini-app
typedef struct gkyl_vlasov_app gkyl_vlasov_app;

/**
 * Construct a new Vlasov mini-app.
 *
 * @param vm App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New vlasov mini-app object.
 */
gkyl_vlasov_app* gkyl_vlasov_app_new(struct gkyl_vm vm);

/**
 * Initialize species and field.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_init_sim(gkyl_vlasov_app* app);

/**
 * Calculate diagnostic moments.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_calc_mom(gkyl_vlasov_app* app);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_mom_write(gkyl_vlasov_app* app, double tm, int frame);

/**
 * Free Vlasov app.
 *
 * @param app App to release.
 */
void gkyl_vlasov_app_release(gkyl_vlasov_app* app);
