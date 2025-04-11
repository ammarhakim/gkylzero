local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.

rho0 = 2.66666666*1.4
rho1 = 1.4
u0 = 1.25
u1 = 0.0
p0 = 4.5
p1 = 1.0

xc = 0.0
yc = 0.0
r = 0.15

-- Simulation parameters.
Nx = 300 -- Cell count (x-direction).
Ny = 300 -- Cell count (y-direction).
Lx = 1.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.25 -- Final simulation time.
num_frames = 5 -- Number of output frames.
field_energy_calcs = INT_MAX -- Number of times to calculate field energy.
integrated_mom_calcs = INT_MAX -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5*Lx, -0.5*Ly },
  upper = { 0.5*Lx, 0.5*Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,
    
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = Euler.new {
      gasGamma = gas_gamma,
      rpType = G0.EulerRP.Roe
    },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local rho = 0.0
      local u = 0.0
      local v = 0.0
      local p = 0.0

      if x < -0.3 then
        rho = rho0 -- Fluid mass density (upper-left).
        u = u0 -- Fluid x-velocity (upper-left).
        p = p0 -- Fluid pressure (upper-left).
      else
        rho = rho1 -- Fluid mass density (lower-left).
        u = u1 -- Fluid x-velocity (lower-left).
        p = p1 -- Fluid pressure (lower-left).
      end

      if ((x - xc)^2 + (y - yc)^2) < r^2 then
        rho = 0.01
        u = 0.0
        p = 0.01
      end
  
      local mom_x = rho*u -- Fluid momentum density (x-direction).
      local mom_y = 0.0 -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = (p/(gas_gamma - 1.0)) + (0.5*rho*u*u) -- Fluid total energy density.
      
      return rho, mom_x, mom_y, mom_z, Etot
    end,
  
    evolve = true -- Evolve species?
  },

  embed_geo = function (t, xn)
    local x, y = xn[1], xn[2]
      
    local phi = 0.0
    if ((x - xc)^2 + (y - yc)^2) < r^2 then
      phi = -1.0
    else
      phi = 1.0
    end

    return phi
  end
}

-- Run application.
momentApp:run()
