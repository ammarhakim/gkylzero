local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left number density.
Tl = 1.0 -- Left temperature.

nr = 0.125 -- Right number density.
Tr = math.sqrt(0.1 / 0.125) -- Right temperature.

vt = 1.0 -- Thermal velocity.
nu = 100.0 -- Collision frequency.

B0 = 1.0 -- Reference magnetic field strength.

-- Simulation parameters.
Nx = 128 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 * vt -- Domain boundary (velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
num_frames = 1 -- Number of output frames.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

pkpmApp = PKPM.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Neutral species.
  neut = PKPM.Species.new {
    modelID = G0.Model.Default,
    charge = charge, mass = mass,
    
    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions (distribution function).
    initDist = function (t, xn)
      local x, vx = xn[1], xn[2]

      local F0 = 0.0
      local T = 0.0
      
      if x < 0.5 then
        F0 = (nl / math.sqrt(2.0 * pi * Tl * Tl)) * (math.exp(-(vx * vx) / (2.0 * Tl * Tl))) -- Distribution function (F0, left).
        T = Tl -- Total temperature (left).
      else
        F0 = (nr / math.sqrt(2.0 * pi * Tr * Tr)) * (math.exp(-(vx * vx) / (2.0 * Tr * Tr))) -- Distribution function (F0right).
        T = Tr -- Total temperature (right).
      end

      local G = (T * T) * F0 -- Distribution function (G).
      
      return F0, G
    end,

    -- Initial conditions (fluid).
    initFluid = function (t, xn)
      local mom_x = 0.0 -- Total momentum density (x-direction).
      local mom_y = 0.0 -- Total momentum density (y-direction).
      local mom_z = 0.0 -- Total momentum density (z-direction).

      return mom_x, mom_y, mom_z
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = B0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()