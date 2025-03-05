local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_pos = 1.0 -- Positron mass.
charge_pos = 1.0 -- Positron charge.

n0 = 1.0 -- Reference number density.
T = 1.0 -- Electron temperature.
Vx_drift = 5.0 -- Drift velocity (x-direction). Normalized to thermal velocity. 

delta_n = 1.0e-4 -- Applied perturbation amplitude.
alphaG = 1.0e-2 -- Normalized ratio of the strength of gravity
                -- compared to the strength of electromagnetism
                -- alphaG = 4 pi G epsilon0 m^2/q^2

-- Derived physical quantities (using normalized code units).
vte = math.sqrt(T / mass_elc) -- Electron thermal velocity.
omega_pe = math.sqrt((charge_elc * charge_elc) * n0 / (epsilon0 * mass_elc)) -- Electron plasma frequency.
lambda_D = vte / omega_pe -- Electron Debye length.
lambda_J = lambda_D/math.sqrt(alphaG) -- Jeans length 
omega_J = vte/lambda_J -- Jeans frequency

-- Lowest perturbed wave number.
-- Needs to be at the largest scale, which is the Jeans length
kx = (0.1 / lambda_J)
mode_init = 1 -- Initial wave mode to perturb with noise.
mode_final = 32 -- Final wave mode to perturb with noise.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 128 -- Cell count (velocity space: vx-direction).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
vx_max = 64.0 * vte -- Domain boundary (velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 / omega_J -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 2,
    projections = {
      -- Two counter-streaming Maxwellians.
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]
          local perturb = 0.0
          math.randomseed(0)
        
          for i = mode_init, mode_final do
            perturb = perturb + delta_n * math.random() * math.cos(i * kx * x + 2.0 * pi * math.random()) 
          end          

          local n = 0.5 * (1.0 + perturb) * n0 -- Total number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift -- Total left-going drift velocity.
        end,

        correctAllMoments = true
      },
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]
          local perturb = 0.0
          math.randomseed(0)
        
          for i = mode_init, mode_final do
            perturb = perturb + delta_n * math.random() * math.cos(i * kx * x + 2.0 * pi * math.random()) 
          end          

          local n = 0.5 * (1.0 + perturb) * n0 -- Total number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return -Vx_drift -- Total right-going drift velocity.
        end,

        correctAllMoments = true
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { "M0", "M1i", "M2" }
  },

  -- Positrons.
  pos = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_pos, mass = mass_pos,
    
    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 2,
    projections = {
      -- Two counter-streaming Maxwellians.
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.5 * n0 -- Total number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift -- Total left-going drift velocity.
        end,

        correctAllMoments = true
      },
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.5 * n0 -- Total number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return -Vx_drift -- Total right-going drift velocity.
        end,

        correctAllMoments = true
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { "M0", "M1i", "M2" }
  },

  isElectrostatic = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0,
    alphaG = alphaG, 

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcPeriodic
      },
      upperType = {
        G0.PoissonBc.bcPeriodic
      }
    }
  }
}

vlasovApp:run()