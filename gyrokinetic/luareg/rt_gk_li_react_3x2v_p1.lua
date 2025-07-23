local Gyrokinetic = G0.Gyrokinetic

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using non-normalized physical units).
epsilon0 = 8.854187817620389850536563031710750260608e-12 -- Permittivity of free space.
mass_elc = 9.10938215e-31 -- Electron mass.
mass_ion = 2.014 * 1.672621637e-27 -- Proton mass.
mass_Li1 = 6.94 * 1.672621637e-27 -- Li1+ ion mass.
mass_Li2 = 6.94 * 1.672621637e-27 -- Li2+ ion mass.
charge_elc = -1.602176487e-19 -- Electron charge.
charge_ion = 1.602176487e-19 -- Proton charge.
charge_Li1 = 1.602176487e-19 -- Li1+ ion charge.
charge_Li2 = 2.0 * 1.602176487e-19 -- Li2+ ion charge.

Te = 100.0 * 1.602176487e-19 -- Electron temperature.
Ti = 100.0 * 1.602176487e-19 -- Ion temperature.
TLi1 = 100.0 * 1.602176487e-19 -- Li1+ temperature.
TLi2 = 100.0 * 1.602176487e-19 -- Li2+ temperature.
n0_elc = 7.0e18 -- Electron reference number density (1 / m^3).

B_axis = 0.5 -- Magnetic field axis (simple toroidal coordinates).
R0 = 0.85 -- Major radius (simple toroidal coordinates).
a0 = 0.15 -- Minor axis (simple toroidal coordinates).

nu_frac = 0.1 -- Collision frequency fraction.

-- Derived physical quantities (using non-normalized physical units).
n0_Li1 = 0.05 * n0_elc -- Li1+ reference number density (1 / m^3).
n0_Li2 = 0.05 * n0_elc -- Li2+ reference number density (1 / m^3).
n0_ion = n0_elc - n0_Li1 - (2.0 * n0_Li2) -- Ion reference number density (1 / m^3).

R = R0 + a0 -- Radial coordinate (simple toroidal coordinates).
B0 = B_axis * (R0 / R) -- Reference magnetic field strength (Tesla).

log_lambda_elc = 6.6 - 0.5 * math.log(n0_elc / 1.0e20) + 1.5 * math.log(Te / charge_ion) -- Electron Coulomb logarithm.
log_lambda_ion = 6.6 - 0.5 * math.log(n0_elc / 1.0e20) + 1.5 * math.log(Ti / charge_ion) -- Ion Coulomb logarithm.
nu_elc = nu_frac * log_lambda_elc * math.pow(charge_ion, 4.0) * n0_elc /
  (6.0 * math.sqrt(2.0) * math.pow(pi, 3.0 / 2.0) * math.pow(epsilon0, 2.0) * math.sqrt(mass_elc) * math.pow(Te, 3.0 / 2.0)) -- Electron collision frequency.
nu_ion = nu_frac * log_lambda_ion * math.pow(charge_ion, 4.0) * n0_elc /
  (12.0 * math.pow(pi, 3.0 / 2.0) * math.pow(epsilon0, 2.0) * math.sqrt(mass_ion) * math.pow(Ti, 3.0 / 2.0)) -- Ion collision frequency.

c_s = math.sqrt(Te / mass_ion) -- Sound speed.
vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
vtLi1 = math.sqrt(TLi1 / mass_Li1) -- Li1+ thermal velocity.
vtLi2 = math.sqrt(TLi2 / mass_Li2) -- Li2+ thermal velocity.
omega_ci = math.abs(charge_ion * B0 / mass_ion) -- Ion cyclotron frequency.
rho_s = c_s / omega_ci -- Ion-sound gyroradius.

n_src = 1.4690539 * 3.612270e23 -- Source number density.
T_src = 2.0 * Te -- Source temperature.
xmu_src = R -- Source mean position (x-direction).
xsigma_src = 0.005 -- Source standard deviation (x-direction).
floor_src = 0.1 -- Minimum source intensity.

-- Simulation parameters.
Nx = 4 -- Cell count (configuration space: x-direction).
Ny = 1 -- Cell count (configuration space: y-direction).
Nz = 8 -- Cell count (configuration space: z-direction).
Nvpar = 6 -- Cell count (velocity space: parallel velocity direction).
Nmu = 4 -- Cell count (velocity space: magnetic moment direction).
Lx = 50.0 * rho_s -- Domain size (configuration space: x-direction).
Ly = 100.0 * rho_s -- Domain size (configuration space: y-direction).
Lz = 4.0 -- Domain size (configuration space: z-direction).
vpar_max_elc = 4.0 * vte -- Domain boundary (electron velocity space: parallel velocity direction).
mu_max_elc = (3.0 / 2.0) * 0.5 * mass_elc * math.pow(4.0 * vte, 2.0) / (2.0 * B0) -- Domain boundary (electron velocity space: magnetic moment direction).
vpar_max_ion = 4.0 * vti -- Domain boundary (ion velocity space: parallel velocity direction).
mu_max_ion = (3.0 / 2.0) * 0.5 * mass_ion * math.pow(4.0 * vti, 2.0) / (2.0 * B0) -- Domain boundary (ion velocity space: magnetic moment direction).
vpar_max_Li1 = 4.0 * vtLi1 -- Domain boundary (Li1+ velocity space: parallel velocity direction).
mu_max_Li1 = (3.0 / 2.0) * 0.5 * mass_ion * math.pow(4.0 * vtLi1, 2.0) / (2.0 * B0) -- Domain boundary (Li1+ velocity space: magnetic moment direction).
vpar_max_Li2 = 4.0 * vtLi2 -- Domain boundary (Li2+ velocity space: parallel velocity direction).
mu_max_Li2 = (3.0 / 2.0) * 0.5 * mass_ion * math.pow(4.0 * vtLi2, 2.0) / (2.0 * B0) -- Domain boundary (Li2+ velocity space: magnetic moment direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1.0e-7 -- Final simulation time.
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
  lower = { R - 0.5 * Lx, -0.5 * Ly, -0.5 * Lz },
  upper = { R + 0.5 * Lx, 0.5 * Ly, 0.5 * Lz },
  cells = { Nx, Ny, Nz },
  cflFrac = cfl_frac,

  --basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions (y-direction only).

  geometry = {
    geometryID = G0.Geometry.MapC2P,
    world = { 0.0, 0.0, 0.0 },

    -- Computational coordinates (x, y, z) from physical coordinates (X, Y, Z).
    mapc2p = function (t, zc)
      local x, y, z = zc[1], zc[2], zc[3]

      local xp = { }
      
      local R = x
      local phi = z / (R0 + a0)
      local X = R * math.cos(phi)
      local Y = R * math.sin(phi)
      local Z = y

      xp[1] = X
      xp[2] = Y
      xp[3] = Z

      return xp[1], xp[2], xp[3]
    end,

    -- Magnetic field strength.
    bmagFunc = function (t, zc)
      local x = zc[1]

      return B0 * R / x
    end
  },

  -- Electrons.
  elc = Gyrokinetic.Species.new {
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vpar_max_elc, 0.0 },
    upper = { vpar_max_elc, mu_max_elc },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0_elc,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        local x, z = xn[1], xn[3]

        local src_density = math.max(math.exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))), floor_src) * n_src
        local src_temp = 0.0
        local n = 0
      
        if x < xmu_src + 3.0 * xsigma_src then
          src_temp = T_src
        else
          src_temp = (3.0 / 8.0) * T_src
        end
      
        local c_s_src = math.sqrt((5.0 / 3.0) * src_temp / mass_ion)
        local n_peak = 4.0 * math.sqrt(5.0) / 3.0 / c_s_src * (0.125 * Lz) * src_density
      
        if math.abs(z) <= 0.25 * Lz then
          n = 0.5 * n_peak * (1.0 + math.sqrt(1.0 - (z / (0.25 * Lz)) * (z / (0.25 * Lz)))) -- Electron total number density (left).
        else
          n = 0.5 * n_peak -- Electron total number density (right).
        end
        
        return n
      end,
      temperatureInit = function (t, xn)
        local x = xn[1]

        local T = 0.0

        if x < xmu_src + 3.0 * xsigma_src then
          T = (5.0 / 4.0) * Te -- Electron isotropic temperature (left).
        else
          T = 0.5 * Te -- Electron isotropic temperature (right).
        end

        return T
      end,
      parallelVelocityInit = function (t, xn)
        return 0.0 -- Electron parallel velocity.
      end
    },

    source = {
      sourceID = G0.Source.Proj,
  
      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.MaxwellianPrimitive,

          densityInit = function (t, xn)
            local x, z = xn[1], xn[3]

            local n = 0.0

            if math.abs(z) < 0.25 * Lz then
              n = math.max(math.exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))),
                floor_src) * n_src -- Electron source total number density (left).
            else
              n = 1.0e-40 * n_src -- Electron source total number density (right).
            end

            return n
          end,
          temperatureInit = function (t, xn)
            local x = xn[1]

            local T = 0.0

            if x < xmu_src + 3.0 * xsigma_src then
              T = T_src -- Electron source isotropic temperature (left).
            else
              T = (3.0 / 8.0) * T_src -- Electron source isotropic temperature (right).
            end

            return T -- Electron source isotropic temperature.
          end,
          parallelVelocityInit = function (t, xn)
            return 0.0 -- Electron source parallel velocity.
          end
        }
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_elc
      end,

      numCrossCollisions = 1,
      collideWith = { "ion" }
    },

    reaction = {
      numReactions = 2,
    
      reactionTypes = {
        {
          reactionID = G0.Reaction.Ionization,
          selfType = G0.Self.Electron,
          ionType = G0.Ion.Lithium,

          electronName = "elc",
          ionName = "Li2",
          donorName = "Li1",
        
          chargeState = 1,
          ionMass = mass_Li2,
          electronMass = mass_elc
        },
        {
          reactionID = G0.Reaction.Recombination,
          selfType = G0.Self.Electron,
          ionType = G0.Ion.Lithium,

          electronName = "elc",
          ionName = "Li2",
          receiverName = "Li1",
        
          chargeState = 1,
          ionMass = mass_Li2,
          electronMass = mass_elc
        }
      },
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcZeroFlux
      },
      upper = {
        type = G0.SpeciesBc.bcZeroFlux
      }
    },
    bcz = {
      lower = {
        type = G0.SpeciesBc.bcGkSheath
      },
      upper = {
        type = G0.SpeciesBc.bcGkSheath
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp }
  },

  -- Ions.
  ion = Gyrokinetic.Species.new {
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vpar_max_ion, 0.0 },
    upper = { vpar_max_ion, mu_max_ion },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0_ion,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        local x, z = xn[1], xn[3]

        local src_density = math.max(math.exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))), floor_src) * n_src
        local src_temp = 0.0
        local n = 0
      
        if x < xmu_src + 3.0 * xsigma_src then
          src_temp = T_src
        else
          src_temp = (3.0 / 8.0) * T_src
        end
      
        local c_s_src = math.sqrt((5.0 / 3.0) * src_temp / mass_ion)
        local n_peak = 0.85 * 4.0 * math.sqrt(5.0) / 3.0 / c_s_src * (0.125 * Lz) * src_density
      
        if math.abs(z) <= 0.25 * Lz then
          n = 0.5 * n_peak * (1.0 + math.sqrt(1.0 - (z / (0.25 * Lz)) * (z / (0.25 * Lz)))) -- Ion total number density (left).
        else
          n = 0.5 * n_peak -- Ion total number density (right).
        end
        
        return n
      end,
      temperatureInit = function (t, xn)
        local x = xn[1]

        local T = 0.0

        if x < xmu_src + 3.0 * xsigma_src then
          T = (5.0 / 4.0) * Ti -- Ion isotropic temperature (left).
        else
          T = 0.5 * Ti -- Ion isotropic temperature (right).
        end

        return T
      end,
      parallelVelocityInit = function (t, xn)
        return 0.0 -- Ion parallel velocity.
      end
    },

    source = {
      sourceID = G0.Source.Proj,
  
      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.MaxwellianPrimitive,

          densityInit = function (t, xn)
            local x, z = xn[1], xn[3]

            local n = 0.0

            if math.abs(z) < 0.25 * Lz then
              n = 0.85 * math.max(math.exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))),
                floor_src) * n_src -- Ion source total number density (left).
            else
              n = 0.85 * 1.0e-40 * n_src -- Ion source total number density (right).
            end

            return n
          end,
          temperatureInit = function (t, xn)
            local x = xn[1]

            local T = 0.0

            if x < xmu_src + 3.0 * xsigma_src then
              T = T_src -- Ion source isotropic temperature (left).
            else
              T = (3.0 / 8.0) * T_src -- Ion source isotropic temperature (right).
            end

            return T -- Ion source isotropic temperature.
          end,
          parallelVelocityInit = function (t, xn)
            return 0.0 -- Ion source parallel velocity.
          end
        }
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ion
      end,

      numCrossCollisions = 1,
      collideWith = { "elc" }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcZeroFlux
      },
      upper = {
        type = G0.SpeciesBc.bcZeroFlux
      }
    },
    bcz = {
      lower = {
        type = G0.SpeciesBc.bcGkSheath
      },
      upper = {
        type = G0.SpeciesBc.bcGkSheath
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp }
  },

  -- Li1+ ions.
  Li1 = Gyrokinetic.Species.new {
    charge = charge_Li1, mass = mass_Li1,
    
    -- Velocity space grid.
    lower = { -vpar_max_Li1, 0.0 },
    upper = { vpar_max_Li1, mu_max_Li1 },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0_Li1,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        local x, z = xn[1], xn[3]

        local src_density = math.max(math.exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))), floor_src) * n_src
        local src_temp = 0.0
        local n = 0
      
        if x < xmu_src + 3.0 * xsigma_src then
          src_temp = T_src
        else
          src_temp = (3.0 / 8.0) * T_src
        end
      
        local c_s_src = math.sqrt((5.0 / 3.0) * src_temp / mass_ion)
        local n_peak = 0.05 * 4.0 * math.sqrt(5.0) / 3.0 / c_s_src * (0.125 * Lz) * src_density
      
        if math.abs(z) <= 0.25 * Lz then
          n = 0.5 * n_peak * (1.0 + math.sqrt(1.0 - (z / (0.25 * Lz)) * (z / (0.25 * Lz)))) -- Li1+ ion total number density (left).
        else
          n = 0.5 * n_peak -- Li1+ ion total number density (right).
        end
        
        return n
      end,
      temperatureInit = function (t, xn)
        local x = xn[1]

        local T = 0.0

        if x < xmu_src + 3.0 * xsigma_src then
          T = (5.0 / 4.0) * Ti -- Li1+ ion isotropic temperature (left).
        else
          T = 0.5 * Ti -- Li1+ ion isotropic temperature (right).
        end

        return T
      end,
      parallelVelocityInit = function (t, xn)
        return 0.0 -- Li1+ ion parallel velocity.
      end
    },

    source = {
      sourceID = G0.Source.Proj,
  
      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.MaxwellianPrimitive,

          densityInit = function (t, xn)
            local x, z = xn[1], xn[3]

            local n = 0.0

            if math.abs(z) < 0.25 * Lz then
              n = 0.05 * math.max(math.exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))),
                floor_src) * n_src -- Li1+ ion source total number density (left).
            else
              n = 0.05 * 1.0e-40 * n_src -- Li1+ ion source total number density (right).
            end

            return n
          end,
          temperatureInit = function (t, xn)
            local x = xn[1]

            local T = 0.0

            if x < xmu_src + 3.0 * xsigma_src then
              T = T_src -- Li1+ ion source isotropic temperature (left).
            else
              T = (3.0 / 8.0) * T_src -- Li1+ ion source isotropic temperature (right).
            end

            return T -- Li1+ ion source isotropic temperature.
          end,
          parallelVelocityInit = function (t, xn)
            return 0.0 -- Li1+ ion source parallel velocity.
          end
        }
      }
    },

    reaction = {
      numReactions = 2,
    
      reactionTypes = {
        {
          reactionID = G0.Reaction.Ionization,
          selfType = G0.Self.Donor,
          ionType = G0.Ion.Lithium,

          electronName = "elc",
          ionName = "Li2",
          donorName = "Li1",
        
          chargeState = 1,
          ionMass = mass_Li2,
          electronMass = mass_elc
        },
        {
          reactionID = G0.Reaction.Recombination,
          selfType = G0.Self.Receiver,
          ionType = G0.Ion.Lithium,

          electronName = "elc",
          ionName = "Li2",
          receiverName = "Li1",
        
          chargeState = 1,
          ionMass = mass_Li2,
          electronMass = mass_elc
        }
      },
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcZeroFlux
      },
      upper = {
        type = G0.SpeciesBc.bcZeroFlux
      }
    },
    bcz = {
      lower = {
        type = G0.SpeciesBc.bcGkSheath
      },
      upper = {
        type = G0.SpeciesBc.bcGkSheath
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp }
  },

  -- Li2+ ions.
  Li2 = Gyrokinetic.Species.new {
    charge = charge_Li2, mass = mass_Li2,
    
    -- Velocity space grid.
    lower = { -vpar_max_Li2, 0.0 },
    upper = { vpar_max_Li2, mu_max_Li2 },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0_Li2,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        local x, z = xn[1], xn[3]

        local src_density = math.max(math.exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))), floor_src) * n_src
        local src_temp = 0.0
        local n = 0
      
        if x < xmu_src + 3.0 * xsigma_src then
          src_temp = T_src
        else
          src_temp = (3.0 / 8.0) * T_src
        end
      
        local c_s_src = math.sqrt((5.0 / 3.0) * src_temp / mass_ion)
        local n_peak = 0.05 * 4.0 * math.sqrt(5.0) / 3.0 / c_s_src * (0.125 * Lz) * src_density
      
        if math.abs(z) <= 0.25 * Lz then
          n = 0.5 * n_peak * (1.0 + math.sqrt(1.0 - (z / (0.25 * Lz)) * (z / (0.25 * Lz)))) -- Li2+ ion total number density (left).
        else
          n = 0.5 * n_peak -- Li2+ ion total number density (right).
        end
        
        return n
      end,
      temperatureInit = function (t, xn)
        local x = xn[1]

        local T = 0.0

        if x < xmu_src + 3.0 * xsigma_src then
          T = (5.0 / 4.0) * Ti -- Li2+ ion isotropic temperature (left).
        else
          T = 0.5 * Ti -- Li2+ ion isotropic temperature (right).
        end

        return T
      end,
      parallelVelocityInit = function (t, xn)
        return 0.0 -- Li2+ ion parallel velocity.
      end
    },

    source = {
      sourceID = G0.Source.Proj,
  
      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.MaxwellianPrimitive,

          densityInit = function (t, xn)
            local x, z = xn[1], xn[3]

            local n = 0.0

            if math.abs(z) < 0.25 * Lz then
              n = 0.05 * math.max(math.exp(-((x - xmu_src) * (x - xmu_src)) / ((2.0 * xsigma_src) * (2.0 * xsigma_src))),
                floor_src) * n_src -- Li2+ ion source total number density (left).
            else
              n = 0.05 * 1.0e-40 * n_src -- Li2+ ion source total number density (right).
            end

            return n
          end,
          temperatureInit = function (t, xn)
            local x = xn[1]

            local T = 0.0

            if x < xmu_src + 3.0 * xsigma_src then
              T = T_src -- Li2+ ion source isotropic temperature (left).
            else
              T = (3.0 / 8.0) * T_src -- Li2+ ion source isotropic temperature (right).
            end

            return T -- Li2+ ion source isotropic temperature.
          end,
          parallelVelocityInit = function (t, xn)
            return 0.0 -- Li2+ ion source parallel velocity.
          end
        }
      }
    },

    reaction = {
      numReactions = 2,
    
      reactionTypes = {
        {
          reactionID = G0.Reaction.Ionization,
          selfType = G0.Self.Ion,
          ionType = G0.Ion.Lithium,

          electronName = "elc",
          ionName = "Li2",
          donorName = "Li1",
        
          chargeState = 1,
          ionMass = mass_Li2,
          electronMass = mass_elc
        },
        {
          reactionID = G0.Reaction.Recombination,
          selfType = G0.Self.Ion,
          ionType = G0.Ion.Lithium,

          electronName = "elc",
          ionName = "Li2",
          receiverName = "Li1",
        
          chargeState = 1,
          ionMass = mass_Li2,
          electronMass = mass_elc
        }
      },
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcZeroFlux
      },
      upper = {
        type = G0.SpeciesBc.bcZeroFlux
      }
    },
    bcz = {
      lower = {
        type = G0.SpeciesBc.bcGkSheath
      },
      upper = {
        type = G0.SpeciesBc.bcGkSheath
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp }
  },

  -- Field.
  field = Gyrokinetic.Field.new {
    femParBc = G0.ParProjBc.None,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcDirichlet,
        G0.PoissonBc.bcPeriodic
      },
      upperType = {
        G0.PoissonBc.bcDirichlet,
        G0.PoissonBc.bcPeriodic
      },
      lowerValue = {
        0.0
      },
      upperValue = {
        0.0
      }
    }
  }
}

gyrokineticApp:run()