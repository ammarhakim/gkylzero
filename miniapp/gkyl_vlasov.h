#pragma once

#ifndef GKYL_MAX_SPECIES
# define GKYL_MAX_SPECIES 2
#endif

// Parameters for a Vlasov species
struct gkyl_vlasov_species {
    char name[64]; // Species name
    double charge, mass; // Charge and mass
    double lower[3], upper[3]; // Lower, upper bounds of velocity-space
    int cells[3]; // Velocity-space cells

    int evolve; // Evolve species? 1-yes, 0-no

    void *ctx; // Context for initial condition init function
    // Pointer to initialization function
    void (*init)(double t, const double *xn, double *fout, void *ctx);
};

// Parameter for EM field
struct gkyl_em_field {
    double epsilon0;
    double mu0;

    int evolve; // Evolve field? 1-yes, 0-no

    void *ctx; // Context for initial condition init function
    // Pointer to initialization function
    void (*init)(double t, const double *xn, double *fout, void *ctx);    
};

// Top-level app parameters
struct gkyl_vm {
    char name[64]; // Name of app: used as output prefix
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
 * @param vm App inputs. See struct docs. All struct params MUST be initialized
 * @return New vlasov mini-app object.
 */
gkyl_vlasov_app* vlasov_app_new(struct gkyl_vm vm);

/**
 * Free Vlasov app.
 *
 * @param app App to release.
 */
void vlasov_app_release(gkyl_vlasov_app* app);


