local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 1836.153 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.

Te_over_Ti = 4.0 -- Electron temperature / ion temperature.
B0 = 1.0 -- Reference magnetic field strength.

-- Derived physical quantities (using normalized code units).
vte = 1.0 -- Electron thermal velocity.
vti = vte / math.sqrt(Te_over_Ti * mass_ion) -- Ion thermal velocity.
cs = vte / math.sqrt(mass_ion) -- Sound speed.

Vx_drift = 2.0 * cs -- Drift velocity (x-direction).

nu_elc = 1.0e-4 -- Electron collision frequency.
nu_ion = 1.0e-4 / math.sqrt(mass_ion) * (Te_over_Ti * math.sqrt(Te_over_Ti)) -- Ion collision frequency.

-- Simulation parameters.
Nx = 128 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
Lx = 256.0 -- Domain size (configuration space: x-direction).
vx_max_elc = 6.0 * vte -- Domain boundary (electron velocity space: vx-direction).
vx_max_ion = 32.0 * vti -- Domain boundary (ion velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 100.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

pkpmApp = PKPM.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  useExplicitSource = true,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Electrons.
  elc = PKPM.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx },

    -- Initial conditions (distribution function).
    initDist = function (t, xn)
      local vx = xn[2]

      local n = (n0 / math.sqrt(2.0 * pi * vte * vte)) * (math.exp(-(vx * vx) / (2.0 * vte * vte))) -- Electron total number density.
      local T_sq_n = (vte * vte) * n -- Electron temperature squared times electron number density.
      
      return n, T_sq_n
    end,

    -- Initial conditions (fluid).
    initFluid = function (t, xn)
      local x = xn[1]
      
      local mom_x = -n0 * mass_elc * Vx_drift * math.tanh(x) -- Total momentum density (x-direction).
      local mom_y = 0.0 -- Total momentum density (y-direction).
      local mom_z = 0.0 -- Total momentum density (z-direction).

      return mom_x, mom_y, mom_z
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_elc -- Electron collision frequency.
      end
    },

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  },

  -- Ions.
  ion = PKPM.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx },

    -- Initial conditions (distribution function).
    initDist = function (t, xn)
      local vx = xn[2]

      local n = (n0 / math.sqrt(2.0 * pi * vti * vti)) * (math.exp(-(vx * vx) / (2.0 * vti * vti))) -- Ion total number density.
      local T_sq_n = (vti * vti) * n -- Ion temperature squared times ion number density.
      
      return n, T_sq_n
    end,

    -- Initial conditions (fluid).
    initFluid = function (t, xn)
      local x = xn[1]

      local mom_x = -n0 * mass_ion * Vx_drift * math.tanh(x) -- Total momentum density (x-direction).
      local mom_y = 0.0 -- Total momentum density (y-direction).
      local mom_z = 0.0 -- Total momentum density (z-direction).

      return mom_x, mom_y, mom_z
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ion -- Ion collision frequency.
      end
    },

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      local Ex = 0.0 -- External electric field (x-direction).
      local Ey = 0.0 -- External electric field (y-direction).
      local Ez = 0.0 -- External electric field (z-direction).

      local Bx = B0 -- External magnetic field (x-direction).
      local By = 0.0 -- External magnetic field (y-direction).
      local Bz = 0.0 -- External magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz
    end,

    evolve = true, -- Evolve field?
    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy }, -- Copy boundary conditions (x-direction).

    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()