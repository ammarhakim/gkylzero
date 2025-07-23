local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

vt = 1.0 -- Thermal velocity.

-- Simulation parameters.
Nr = 16 -- Cell count (configuration space: radial direction).
Ntheta = 32 -- Cell count (configuration space: angular direction).
Nvr = 8 -- Cell count (velocity space: radial direction).
Nvtheta = 32 -- Cell count (velocity space: angular direction).
Lr = 1.5 -- Domain size (configuration space: radial direction).
Ltheta = 2.0 * pi -- Domain size (configuration space: angular direction).
vr_max = 0.5 * vt -- Domain boundary (velocity space: radial direction).
vtheta_max = 3.0 * vt -- Domain boundary (velocity space: angular direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.01 -- Final simulation time.
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
  lower = { 0.5, 0.0 },
  upper = { 0.5 + Lr, Ltheta },
  cells = { Nr, Ntheta },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions (angular direction only).

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,

    hamiltonian = function (t, xn)
      local q_r, p_r_dot, p_theta_dot = xn[1], xn[3], xn[4]

      local inv_metric_r_r = 1.0
      local inv_metric_r_theta = 0.0
      local inv_metric_theta_theta = 1.0 / (q_r * q_r)

      local hamiltonian = (0.5 * inv_metric_r_r * p_r_dot * p_r_dot) + (0.5 * (2.0 * inv_metric_r_theta * p_r_dot * p_theta_dot)) +
        (0.5 * inv_metric_theta_theta * p_theta_dot * p_theta_dot) - 1.0/q_r -- Canonical Hamiltonian.

      return hamiltonian
    end,
    inverseMetric = function (t, xn)
      local q_r = xn[1]

      local inv_metric_r_r = 1.0 -- Inverse metric tensor (radial-radial component).
      local inv_metric_r_theta = 0.0 -- Inverse metric tensor (radial-angular component).
      local inv_metric_theta_theta = 1.0 / (q_r * q_r) -- Inverse metric tensor (angular-angular component).

      return inv_metric_r_r, inv_metric_r_theta, inv_metric_theta_theta
    end,
    metric = function (t, xn)
      local q_r = xn[1]

      local metric_r_r = 1.0 -- Metric tensor (radial-radial component).
      local metric_r_theta = 0.0 -- Metric tensor (radial-angular component).
      local metric_theta_theta = q_r * q_r -- Metric tensor (angular-angular component).

      return metric_r_r, metric_r_theta, metric_theta_theta
    end,
    metricDeterminant = function (t, xn)
      local q_r = xn[1]

      local metric_det = q_r -- Metric tensor determinant.

      return metric_det
    end,
    
    -- Velocity space grid.
    lower = { -vr_max, 0.5 },
    upper = { vr_max, 1.75 },
    cells = { Nvr, Nvtheta },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)

          local q_r, q_theta, p_r_dot, p_theta_dot = xn[1], xn[2], xn[3], xn[4]
          local metric_det = q_r -- Metric tensor determinant.
          local f = 0.0 -- background distribution function
         
          -- Initalize f's positisional amplitude in space  
          -- Init with a sqaure function in momentum space
          -- around orbital speed p_theta = r (for circular orbits)
          -- and r_dot = 0
          local sigma = 0.15
          if ((p_r_dot > - 0.1) and (p_r_dot < 0.1)) then
            if ((p_theta_dot > math.sqrt(q_r - 0.1)) and (p_theta_dot < math.sqrt(q_r + 0.1))) then
              if (q_r < 1.0) then
                  f = 0.1 + 0.9 * math.exp(-((q_r - 1.0) * (q_r - 1.0)) / (2.0 * sigma * sigma))
              elseif (q_r > 1.5) then
                  f =  0.1 + 0.9 * math.exp(-((q_r - 1.5) * (q_r - 1.5)) / (2.0 * sigma * sigma))
              else
                  f = 1.0
              end
            end
          end
        
          -- Set distribution function, multiply by jacobian for the correct density
          local jac = q_r
          return jac*f
        end
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcAbsorb
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1i_from_H }
  },

  skipField = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0,

    isStatic = true
  }
}

vlasovApp:run()