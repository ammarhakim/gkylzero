local Moments = G0.Moments
local ReactiveEuler = G0.Moments.Eq.ReactiveEuler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.
specific_heat_capacity = 2.5 -- Specific heat capacity.
energy_of_formation = 1.0 -- Energy of formation.
ignition_temperature = 0.25 -- Ignition temperature.
reaction_rate = 250.0 -- Reaction rate.

rhol = 1.4 -- Left fluid mass density.
ul = 0.0 -- Left fluid velocity.
pl = 1.0 -- Left fluid pressure.

rhor = 0.8887565 -- Right fluid mass density.
ur = -0.577350 -- Right fluid velocity.
pr = 0.191709 -- Right fluid pressure.

-- Simulation parameters.
Nx = 512 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.5 -- Final simulation time.
num_frames = 1 -- Number of output frames.

momentApp = Moments.App.new {
  tEnd = t_end,
  nFrame = num_frames,
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,
  
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Fluid.
  fluid = Moments.Species.new {
    equation = ReactiveEuler.new {
      gasGamma = gas_gamma,
      specificHeatCapacity = specific_heat_capacity,
      energyOfFormation = energy_of_formation,
      ignitionTemperature = ignition_temperature,
      reactionRate = reaction_rate
    },

    hasReactivity = true,
    reactivityGasGamma = gas_gamma,
    reactivitySpecificHeatCapacity = specific_heat_capacity,
    reactivityEnergyOfFormation = energy_of_formation,
    reactivityIgnitionTemperature = ignition_temperature,
    reactivityReactionRate = reaction_rate,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local rho = 0.0
      local u = 0.0
      local p = 0.0
  
      if x < 0.25 then
        rho = rhol -- Fluid mass density (left).
        u = ul -- Fluid velocity (left).
        p = pl -- Fluid pressure (left).
      else
        rho = rhor -- Fluid mass density (right).
        u = ur -- Fluid velocity (right).
        p = pr -- Fluid pressure (right).
      end

      local mom_x = rho * u -- Fluid momentum density (x-direction).
      local mom_y = 0.0 -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = (p / (gas_gamma - 1.0)) + (0.5 * rho * u * u) -- Fluid total energy density.
      local reac = rho -- Fluid reaction progress.
    
      return rho, mom_x, mom_y, mom_z, Etot, reac
    end,

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
  },
  
  -- Field.
  field = Moments.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy }, -- Copy boundary conditions (x-direction).
  },
}

-- Run application.
momentApp:run()
