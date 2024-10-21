local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Reference density.
vt = 0.2 -- Thermal velocity.
Vx_drift = 1.0 -- Drift velocity (x-direction).

alpha = 1.0e-6 -- Applied perturbation amplitude.
kx = 0.5 -- Perturbed wave number (x-direction).

-- Derived physical quantities (using normalized code units).
T = (vt * vt) * mass_elc -- Temperature.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: vx-direction).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
vx_max = 6.0 -- Domain boundary (velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.6 -- CFL coefficient.

t_end = 40.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
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
        projectionID = "LTE",

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.5 * (1.0 + alpha * math.cos(kx * x)) * n0 -- Total number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T -- Total temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift -- Total left-going drift velocity.
        end,

        correctAllMoments = true
      },
      {
        projectionID = "LTE",

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.5 * (1.0 + alpha * math.cos(kx * x)) * n0 -- Total number density.
          return n
        end,
        temperatureInit = function (t, xn)
          return T -- Total temperature.
        end,
        driftVelocityInit = function (t, xn)
          return -Vx_drift -- Total right-going drift velocity.
        end,

        correctAllMoments = true
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { "M0", "M1i", "M2"  }
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = -alpha * math.sin(kx * x) / kx -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,
    
    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()