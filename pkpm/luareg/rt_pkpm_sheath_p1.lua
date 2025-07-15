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

Te_over_Ti = 1.0 -- Electron temperature / ion temperature.
B0 = 1.0 -- Reference magnetic field strength.

-- Derived physical quantities (using normalized code units).
vte = 1.0 -- Electron thermal velocity.
vti = vte / math.sqrt(Te_over_Ti * mass_ion) -- Ion thermal velocity.

nu_elc = 1.0e-5 -- Electron collision frequency.
nu_ion = 1.0e-5 / math.sqrt(mass_ion) * (Te_over_Ti * math.sqrt(Te_over_Ti)) -- Ion collision frequency.

-- Simulation parameters.
Nx = 128 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: vx-direction).
Lx = 256.0 -- Domain size (configuration space: x-direction).
vx_max_elc = 8.0 * vte -- Domain boundary (electron velocity space: vx-direction).
vx_max_ion = 6.0 * vti -- Domain boundary (ion velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 100.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Training parameters.
train_nn = false -- Train neural network on simulation data?
train_ab_initio = true -- Train neural network ab initio?
nn_width = 256 -- Number of neurons to use per layer.
nn_depth = 5 -- Number of layers to use.
train_nn_file = "rt_pkpm_sheath_p1_moms_nn_1" -- File path of neural network to train.
num_trains = GKYL_MAX_INT -- Number of times to train neural network.
num_nn_writes = 1 -- Number of times to write out neural network.
input_moms = { 1, 3, 4 } -- Array of "input" moments to train on.
output_moms = { 5, 6 } -- Array of "output" moments to train on.
test_nn = false -- Test neural network on simulation data?
test_nn_file = "rt_pkpm_sheath_p1_moms_nn_1" -- File path of neural network to test.
num_tests = 1 -- Number of times to test neural network.

pkpmApp = PKPM.App.new {

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

      local F0 = (n0 / math.sqrt(2.0 * pi * vte * vte)) * (math.exp(-(vx * vx) / (2.0 * vte * vte))) -- Electron distribution function (F0).
      local G = (vte * vte) * F0 -- Electron distribution function (G).
      
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
        return nu_elc -- Electron collision frequency.
      end
    },

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcAbsorb, G0.SpeciesBc.bcAbsorb } -- Absorbing boundary conditions (x-direction).
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

      local F0 = (n0 / math.sqrt(2.0 * pi * vti * vti)) * (math.exp(-(vx * vx) / (2.0 * vti * vti))) -- Ion distribution function (F0).
      local G = (vti * vti) * F0 -- Ion distribution function (G).
      
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
        return nu_ion -- Ion collision frequency.
      end
    },

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcAbsorb, G0.SpeciesBc.bcAbsorb } -- Absorbing boundary conditions (x-direction).
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
    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall }, -- PEC wall boundary conditions (x-direction).
    
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  },

  -- Training parameters.
  trainNN = train_nn,
  trainAbInitio = train_ab_initio,
  NNWidth = nn_width,
  NNDepth = nn_depth,
  trainNNFile = train_nn_file,
  numTrains = num_trains,
  numNNWrites = num_nn_writes,
  inputMoms = input_moms,
  outputMoms = output_moms,
  testNN = test_nn,
  testNNFile = test_nn_file,
  numTests = num_tests
}

pkpmApp:run()