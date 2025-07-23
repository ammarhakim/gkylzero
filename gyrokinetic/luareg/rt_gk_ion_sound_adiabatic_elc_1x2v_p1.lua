local Gyrokinetic = G0.Gyrokinetic

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass_ion = 1.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.
mass_elc = 1.0 / 1836.16 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

Te = 1.0 -- Electron temperature.
Ti = 1.0 -- Ion temperature.
n0 = 1.0 -- Reference number density.
B0 = 1.0 -- Reference magnetic field strength.

k_perp = 0.0 -- Perpendicular wave number (for Poisson solver).

alpha = 0.01 -- Applied perturbation amplitude.
kz = 0.5 -- Perturbed wave number (z-direction).

-- Derived physical quantities (using non-normalized physical units).
nu_ion = 0.005 -- Ion collision frequency.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.

-- Simulation parameters.
Nz = 32 -- Cell count (configuration space: z-direction).
Nvpar = 48 -- Cell count (velocity space: parallel velocity direction).
Nmu = 8 -- Cell count (velocity space: magnetic moment direction).
Lz = 2.0 * pi / kz -- Domain size (configuration space: z-direction).
vpar_max_ion = 6.0 * vti -- Domain boundary (ion velocity space: parallel velocity direction).
mu_max_ion = mass_ion * math.pow(5.0 * vti, 2.0) / (2.0 * B0) -- Domain boundary (ion velocity space: magnetic moment direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

gyrokineticApp = Gyrokinetic.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lz },
  upper = { 0.5 * Lz },
  cells = { Nz },
  cflFrac = cfl_frac,

  --basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  geometry = {
    geometryID = G0.Geometry.MapC2P,
    world = { 0.0, 0.0, 0.0 },

    -- Computational coordinates (x, y, z) from physical coordinates (X, Y, Z).
    mapc2p = function (t, zc)
      local xp = { }
      
      xp[1] = zc[1]
      xp[2] = zc[2]
      xp[3] = zc[3]

      return xp[1], xp[2], xp[3]
    end,

    -- Magnetic field strength.
    bmagFunc = function (t, zc)
      return B0
    end
  },

  -- Ions.
  ion = Gyrokinetic.Species.new {
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vpar_max_ion, 0.0 },
    upper = { vpar_max_ion, mu_max_ion },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        local z = xn[1]

        local n = (1.0 + alpha * math.cos(kz * z)) * n0 -- Ion total number density.

        return n
      end,
      temperatureInit = function (t, xn)
        return Ti -- Ion isotropic temperature.
      end,
      parallelVelocityInit = function (t, xn)
        return 0.0 -- Ion parallel velocity.
      end
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ion
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp }
  },

  -- Field.
  field = Gyrokinetic.Field.new {
    fieldID = G0.GKField.Adiabatic,

    electronMass = mass_elc,
    electronCharge = charge_elc,
    electronDensity = n0,
    electronTemperature = Te,

    femParBc = G0.ParProjBc.Periodic,
    kPerpSq = k_perp * k_perp
  }
}

gyrokineticApp:run()