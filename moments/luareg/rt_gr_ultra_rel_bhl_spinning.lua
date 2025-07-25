-- 2D Bondi-Hoyle-Lyttleton ultra-relativistic accretion problem onto a non-static (Kerr) black hole, for the general relativistic Euler equations,
-- assuming a stiff equation of state.
-- Input parameters describe wind accretion of an ultra-relativistic gas onto a spinning black hole.
-- Based on the analytical solution for stiff relativistic fluids presented in the article:
-- L. I. Petrich, S. L. Shapiro and S. A. Teukolsky (1988), "Accretion onto a moving black hole: An exact solution",
-- Physical Review Letters, Volume 60 (18): 1781-1784.
-- https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.60.1781

local Moments = G0.Moments
local GRUltraRelEuler = G0.Moments.Eq.GRUltraRelEuler
local BlackHole = G0.Moments.Spacetime.BlackHole

-- Physical constants (using normalized code units).
gas_gamma = 2.0 -- Adiabatic index.

rhol = 10.0 -- Left fluid mass density.
ul = 0.3 -- Left fluid velocity.

rhor = 10.0 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.

-- Spacetime parameters (using geometric units).
mass = 0.3 -- Mass of the black hole.
spin = -0.6 -- Spin of the black hole.

pos_x = 2.5 -- Position of the black hole (x-direction).
pos_y = 2.5 -- Position of the black hole (y-direction).
pos_z = 0.0 -- Position of the black hole (z-direction).

-- Simulation parameters.
Nx = 256 -- Cell count (x-direction).
Ny = 256 -- Cell count (y-direction).
Lx = 5.0 -- Domain size (x-direction).
Ly = 5.0 -- Domain size (y-direction).
cfl_frac = 0.95 -- CFL coefficient.

reinit_freq = 100 -- Spacetime reinitialization frequency.

t_end = 15.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

x_loc = 1.0 -- Shock location (x-direction).

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Fluid.
  fluid = Moments.Species.new {
    equation = GRUltraRelEuler.new {
      gasGamma = gas_gamma,
      blackHoleParameters = {
        mass = mass,
        spin = spin,
        posX = pos_x,
        posY = pos_y,
        posZ = pos_z
      },
      reinitFreq = reinit_freq
    },

    hasGRUltraRel = true,
    GRUltraRelGasGamma = gas_gamma,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]
      
      local rho = 0.0
      local u = 0.0

      if x < x_loc then
        rho = rhol -- Fluid mass density (left).
        u = ul -- Fluid velocity (left).
      else
        rho = rhor -- Fluid mass density (right).
        u = ur -- Fluid velocity (right).
      end

      local lapse = BlackHole.lapseFunction(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local shift = BlackHole.shiftVector(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local spatial_metric = BlackHole.spatialMetricTensor(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local spatial_det = BlackHole.spatialMetricDeterminant(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local extrinsic_curvature = BlackHole.extrinsicCurvatureTensor(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local in_excision_region = BlackHole.excisionRegion(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      
      local lapse_der = BlackHole.lapseFunctionDer(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local shift_der = BlackHole.shiftVectorDer(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local spatial_metric_der = BlackHole.spatialMetricTensorDer(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))

      local vel = { u, 0.0, 0.0 }
      local v_sq = 0.0

      for i = 1, 3 do
        for j = 1, 3 do
          v_sq = v_sq + spatial_metric[i][j] * vel[i] * vel[j]
        end
      end

      local W = 1.0 / math.sqrt(1.0 - v_sq)
      if v_sq > 1.0 - math.pow(10.0, -8.0) then
        W = 1.0 / math.sqrt(1.0 - pow(10.0, -8.0))
      end

      local p = (gas_gamma - 1.0) * rho

      local Etot = math.sqrt(spatial_det) * (((rho + p) * (W * W)) - p) -- Fluid total energy density.
      local mom_x = math.sqrt(spatial_det) * (rho + p) * (W * W) * u -- Fluid momentum density (x-direction).
      local mom_y = 0.0 -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).

      local excision = 0.0
      if in_excision_region then
        Etot, mom_x, mom_y, mom_z, lapse = 0.0, 0.0, 0.0, 0.0, 0.0
        for i = 1, 3 do
          shift[i] = 0.0
          lapse_der[i] = 0.0
          for j = 1, 3 do
            spatial_metric[i][j] = 0.0
            extrinsic_curvature[i][j] = 0.0
            shift_der[i][j] = 0.0
            for k = 1, 3 do
              spatial_metric_der[i][j][k] = 0.0
            end
          end  
        end
        
        excision = -1.0
      else
        excision = 1.0
      end
    
      return Etot, mom_x, mom_y, mom_z,
        lapse,
        shift[1], shift[2], shift[3],
        spatial_metric[1][1], spatial_metric[1][2], spatial_metric[1][3],
        spatial_metric[2][1], spatial_metric[2][2], spatial_metric[2][3],
        spatial_metric[3][1], spatial_metric[3][2], spatial_metric[3][3],
        extrinsic_curvature[1][1], extrinsic_curvature[1][2], extrinsic_curvature[1][3],
        extrinsic_curvature[2][1], extrinsic_curvature[2][2], extrinsic_curvature[2][3],
        extrinsic_curvature[3][1], extrinsic_curvature[3][2], extrinsic_curvature[3][3],
        excision,
        lapse_der[1], lapse_der[2], lapse_der[3],
        shift_der[1][1], shift_der[1][2], shift_der[1][3],
        shift_der[2][1], shift_der[2][2], shift_der[2][3],
        shift_der[3][1], shift_der[3][2], shift_der[3][3],
        spatial_metric_der[1][1][1], spatial_metric_der[1][1][2], spatial_metric_der[1][1][3],
        spatial_metric_der[1][2][1], spatial_metric_der[1][2][2], spatial_metric_der[1][2][3],
        spatial_metric_der[1][3][1], spatial_metric_der[1][3][2], spatial_metric_der[1][3][3],
        spatial_metric_der[2][1][1], spatial_metric_der[2][1][2], spatial_metric_der[2][1][3],
        spatial_metric_der[2][2][1], spatial_metric_der[2][2][2], spatial_metric_der[2][2][3],
        spatial_metric_der[2][3][1], spatial_metric_der[2][3][2], spatial_metric_der[2][3][3],
        spatial_metric_der[3][1][1], spatial_metric_der[3][1][2], spatial_metric_der[3][1][3],
        spatial_metric_der[3][2][1], spatial_metric_der[3][2][2], spatial_metric_der[3][2][3],
        spatial_metric_der[3][3][1], spatial_metric_der[3][3][2], spatial_metric_der[3][3][3],
        0.0,
        x, y, 0.0
    end,

    evolve = true, -- Evolve species?
    forceLowOrderFlux = true, -- Use Lax fluxes.
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
