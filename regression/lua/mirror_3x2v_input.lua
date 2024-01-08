-- Maxwell Rosen 10/24/2023
-- Working 3x2v gyrokinetic calculation in a mirror geometry modeled after
-- WHAM.
-- Coordinates are psi,theta,z to be field aligned. May need modification for using psi
--
-- Big changes
-- I changed lowercase z to lineLength
-- Boundary condition bcx and bcz in the species
-- Initial conditions will be updated with loading from file
-- BC on psi? I need periodic in angle?
local Plasma    = require("App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local Logger    = require "Lib.Logger"
-- SciLua (and extensions) used for AD, root finding and integration.
-- Use sci.math instead of built-in 'math', see https://scilua.org/sci_diff.html.
local math      = require("sci.math").generic
local root      = require "sci.root"
local quad      = require "sci.quad"
local diff      = require "sci.diff-recursive"
local df        = diff.df
local log       = Logger { logToFile = true }
log("\n")
log("  3x2v gyrokinetic WHAM simulation\n")
log("\n")
polyOrder    = 1
-- Universal constant parameters.
eps0, mu0    = Constants.EPSILON0, Constants.MU0
eV           = Constants.ELEMENTARY_CHARGE
qe, qi       = -eV, eV
me, mp       = Constants.ELECTRON_MASS, Constants.PROTON_MASS
-- Plasma parameters.
mi           = 2.014 * mp -- Deuterium ion mass.
Te0          = 940 * eV
n0           = 3.e19
B_p          = 0.53
beta         = 0.4                                         -- Ratio of plasma to magnetic pressure.
tau          = (B_p ^ 2) * beta / (2 * mu0 * n0 * Te0) - 1 -- Ti/Te ratio.
Ti0          = tau * Te0
kperpRhos    = 0.1                                         -- k_perp*rho_s in the Poisson equation.
-- Parameters controlling initial conditions.
alim         = 0.125 -- Plasma limiting radius.
alphaIC      = { 2, 10 }
nuFrac       = 1.0 -- Factor multiplying the collisionality.
   else
-- Electron-electron collision freq.
logLambdaElc = 6.6 - 0.5 * math.log(n0 / 1e20) + 1.5 * math.log(Te0 / eV)
nuElc        = nuFrac * logLambdaElc * (eV ^ 4) * n0 /
    (6 * math.sqrt(2) * (math.pi ^ (3 / 2)) * (eps0 ^ 2) * math.sqrt(me) * (Te0 ^ (3 / 2)))
-- Ion-ion collision freq.
logLambdaIon = 6.6 - 0.5 * math.log(n0 / 1e20) + 1.5 * math.log(Ti0 / eV)
nuIon        = nuFrac * logLambdaIon * (eV ^ 4) * n0 /
    (12 * (math.pi ^ (3 / 2)) * (eps0 ^ 2) * math.sqrt(mi) * (Ti0 ^ (3 / 2)))

log(string.format("  Collision frequencies (including nuFrac=%f):\n", nuFrac))
log(string.format("    nu_ee = %e Hz\n", nuElc))
log(string.format("    nu_ii = %e Hz\n", nuIon))
log("\n")

-- Thermal speeds.
vti                = math.sqrt(Ti0 / mi)
vte                = math.sqrt(Te0 / me)
c_s                = math.sqrt(Te0 / mi)

-- Gyrofrequencies and gyroradii.
omega_ci           = eV * B_p / mi
rho_s              = c_s / omega_ci

-- Geometry parameters.
numCellLineLength  = 36 --288
numCellPsi         = 15 -- 36
numCellTheta       = 15 -- 36
-- Axial coordinate Z extents. Endure that Z=0 is not on
-- the boundary of a cell (due to AD errors).
lowerZ      = -2.50
upperZ      = 2.50
RminAtZ0, RmaxAtZ0 = 0.01, 0.10 -- Radial extents

-- Parameters controlling the magnetic equilibrium model.
eqModel            = {
   mcB   = 6.51292,
   gamma = 0.124904,
   Z_m   = 0.98,
}

local function psi_RZ(RIn, ZIn)
   -- Compute flux function value psi as a function of R and Z.
   -- The input eqModel is a dictionary of model parameters.
   local mcB   = eqModel["mcB"]
   local gamma = eqModel["gamma"]
   local Z_m   = eqModel["Z_m"]

   local psi   = 0.5 * (RIn ^ 2) * mcB * (1. / (math.pi * gamma * (1. + ((ZIn - Z_m) / gamma) ^ 2))
      + 1. / (math.pi * gamma * (1. + ((ZIn + Z_m) / gamma) ^ 2)))
   return psi
end

-- Compute the radius R as a function of psi and Z.
-- The input eqModel is a dictionary of model parameters.
local function R_psiZ(psiIn, ZIn)
   local mcB   = eqModel["mcB"]
   local gamma = eqModel["gamma"]
   local Z_m   = eqModel["Z_m"]

   local Rout  = math.sqrt(2. * psiIn / (mcB * (1. / (math.pi * gamma * (1. + ((ZIn - Z_m) / gamma) ^ 2))
      + 1. / (math.pi * gamma * (1. + ((ZIn + Z_m) / gamma) ^ 2)))))
   return Rout
end

local function Bfield_psiZ(psiIn, ZIn)
   local Rcoord = R_psiZ(psiIn, ZIn)

   local mcB    = eqModel["mcB"]
   local gamma  = eqModel["gamma"]
   local Z_m    = eqModel["Z_m"]

   local BRad   = -(1. / 2.) * Rcoord * mcB *
       (-2. * (ZIn - Z_m) / (math.pi * (gamma ^ 3) * ((1. + ((ZIn - Z_m) / gamma) ^ 2) ^ 2))
          - 2. * (ZIn + Z_m) / (math.pi * (gamma ^ 3) * ((1. + ((ZIn + Z_m) / gamma) ^ 2) ^ 2)))

   local BZ     = mcB * (1. / (math.pi * gamma * (1. + ((ZIn - Z_m) / gamma) ^ 2))
      + 1. / (math.pi * gamma * (1. + ((ZIn + Z_m) / gamma) ^ 2)))

   local Bmag   = math.sqrt(BRad ^ 2 + BZ ^ 2)

   return BRad, BZ, Bmag
end

-- Coordinate along the field line (lineLength) as a function of the axial coordinate (Z).
local function lineLength_psiZ(psiIn, ZIn)
   local function integrand(t)
      local _, B_Z, Bmag = Bfield_psiZ(psiIn, t)
      return Bmag / B_Z
   end
   local integral
   local eps = 0.0
   if diff.le(eps, ZIn) then
      integral, _ = quad.dblexp(integrand, eps, ZIn, 1e-14)
   else
      integral, _ = -quad.dblexp(integrand, ZIn, eps, 1e-14)
   end
   return integral
end


-- Stopping criteria for root-finding.
local function stop(tol)
   return function(x, y, xl, xu, yl, yu)
      if diff.lt(math.abs(y), tol) then
         return true
      else
         return false
      end
   end
end

-- Invert lineLength(Z) via root-finding.
local function Z_psiLineLength(psiIn, lineLengthIn)
   local function lossF(Z)
      -- Function we will find the roots off. The root corresponds to
      -- the value of Z that satisfies our coordinate transformation
      return lineLengthIn - lineLength_psiZ(psiIn, Z)
   end
   local maxL = upperZ - lowerZ
   local Zout = 0.0
   local eps  = maxL /
       numCellLineLength -- Interestingly using a smaller eps yields larger errors in some geo quantities.
   if diff.le(0., lineLengthIn) then
      Zout = root.ridders(lossF, -eps, upperZ + maxL / numCellLineLength, stop(1e-14))
   else
      Zout = root.ridders(lossF, lowerZ - maxL / numCellLineLength, eps, stop(1e-14))
   end
   return Zout
end

psiMin = psi_RZ(RminAtZ0, 0)
psiMax = psi_RZ(RmaxAtZ0, 0)

lineLengthMin = lineLength_psiZ(psiMax, lowerZ)
lineLengthMax = lineLength_psiZ(psiMax, upperZ)

local nl                     = 1024
-- Find the throat location with a brute-force search using a finer mesh.
-- Used for initial conditions. Not generalized to 3D.
local B_m, Z_m, lineLength_m = -1.e19, 0.0, 0.0
local psiAvg                 = (psiMax - psiMin) / 2
for k = 1, nl + 1 do
   local L          = lineLengthMax - lineLengthMin       -- Size of the box along the field line.
   local lineLength = lineLengthMin + (k - 1) * (L / nl)  -- Field-aligned coordinate.
   local Z          = Z_psiLineLength(psiAvg, lineLength) -- Cylindrical axial coordinate Z.
   local _, _, Bmag = Bfield_psiZ(psiAvg, Z)
   if B_m < Bmag then
      B_m = Bmag
      Z_m = math.abs(Z)
      lineLength_m = math.abs(lineLength)
   end
end

-- Locate the banana tip, assuming particles injected at
-- 45 degree angle at Z=0 (and that all particles have the same
-- pitch angle) so the banatip occurs where B/B(Z=0)=2.
local BmagDiff = 1.e19
local R_bt, Z_bt, lineLength_bt, B_bt = -1.e19, -1.e19, -1.e19, -1.e19
local _, _, BmagAtReq0 = Bfield_psiZ(psiAvg, 0.0)
for k = 1, nl + 1 do
   local Lp         = 2 * lineLength_m                    -- Size of the plasma along the field line.
   local lineLength = -lineLength_m + (k - 1) * (Lp / nl) -- Field-aligned coordinate.
   local Z          = Z_psiLineLength(psiAvg, lineLength) -- Cylindrical axial coordinate Z.
   local R          = R_psiZ(psiAvg, Z)                   -- Cylindrical radial coordinate R.
   local _, _, Bmag = Bfield_psiZ(psiAvg, Z)
   local currBdiff  = math.abs(Bmag - 2. * BmagAtReq0)
   if currBdiff < BmagDiff then
      BmagDiff = currBdiff
      R_bt = R
      Z_bt = math.abs(Z)
      lineLength_bt = math.abs(lineLength)
      B_bt = Bmag
   end
end

-- Print banana tip location and magnetic field there.
log("  Banana tip field & location:\n")
log(string.format("    B_bt = %f T\n", B_bt))
log(string.format("    R_bt = %f m\n", R_bt))
log(string.format("    Z_bt = %f m\n", Z_bt))
log(string.format("    lineLength_bt = %f m\n", lineLength_bt))

-- Some physics parameters at mirror throats.
local n_m, Te_m = 0.0, 0.0
for k = 1, nl + 1 do
   local L          = lineLengthMax - lineLengthMin       -- Size of the box along the field line.
   local lineLength = lineLengthMin + (k - 1) * (L / nl)  -- Field-aligned coordinate.
   local Z          = Z_psiLineLength(psiAvg, lineLength) -- Cylindrical axial coordinate Z.
   local R          = R_psiZ(psiAvg, Z)                   -- Cylindrical radial coordinate R.
   local _, _, Bmag = Bfield_psiZ(psiAvg, Z)
   if math.abs(B_m - Bmag) < 1.e-10 then
      n_m  = n0 * math.sqrt((1. - ((R - R_bt) / alim) ^ 2) ^ alphaIC[2])
      Te_m = Te0 * math.sqrt((1. - ((R - R_bt) / alim) ^ 2) ^ alphaIC[2])
   end
end
local Ti_m = tau * Te_m
local cs_m = math.sqrt(Te_m * (1. + tau) / mi)

-- Print mirror throat location and magnetic field there.
log("  At mirror throat:\n")
log(string.format("    B_m  = %f T\n", B_m))
log(string.format("    Z_m  = %f m\n", Z_m))
log(string.format("    lineLength_m  = %f m\n", lineLength_m))
log(string.format("    n_m  = %e m^{-3}\n", n_m))
log(string.format("    Te_m = %f eV\n", Te_m / eV))
log(string.format("    Ti_m = %f eV\n", Ti_m / eV))
log(string.format("    cs_m = %e m/s\n", cs_m))

-- Source parameters.
-- Should these all be local? Debugger flag is going off.
--These are global variables that are lowercase
NSrcIon              = 3.1715e23 / 8.
lineLengthSrcIon     = 0.0
sigSrcIon            = Z_m / 4.
NSrcFloorIon         = 0.05 * NSrcIon
TSrc0Ion             = Ti0 * 1.25
TSrcFloorIon         = TSrc0Ion / 8.
NSrcElc              = NSrcIon
lineLengthSrcElc     = lineLengthSrcIon
sigSrcElc            = sigSrcIon
NSrcFloorElc         = NSrcFloorIon
TSrc0Elc             = TSrc0Ion / tau
TSrcFloorElc         = TSrcFloorIon / tau

-- Source functions.
srcDenElc            = function(t, xn)
   local psi, lineLength     = xn[1], xn[3]
   local Z                   = Z_psiLineLength(psi, lineLength) -- Cylindrical axial coordinate.
   local NSrc, lineLengthSrc = NSrcElc, lineLengthSrcElc
   local sigSrc, NSrcFloor   = sigSrcElc, NSrcFloorElc
   if math.abs(Z) <= Z_m then
      return math.max(NSrcFloor, (NSrc / math.sqrt(2. * math.pi * (sigSrc ^ 2))) *
         math.exp(-((lineLength - lineLengthSrc) ^ 2) / (2. * (sigSrc ^ 2))))
   else
      return 1.e-16
   end
end
srcTempElc           = function(t, xn)
   local lineLength = xn[3]
   local sigSrc = sigSrcElc
   local TSrc0, Tfloor = TSrc0Elc, TSrcFloorElc
   if math.abs(lineLength) <= 2. * sigSrc then
      return TSrc0
   else
      return Tfloor
   end
end
srcDenIon            = function(t, xn)
   local psi, lineLength     = xn[1], xn[3]
   local Z                   = Z_psiLineLength(psi, lineLength) -- Cylindrical axial coordinate.
   local NSrc, lineLengthSrc = NSrcIon, lineLengthSrcIon
   local sigSrc, NSrcFloor   = sigSrcIon, NSrcFloorIon
   if math.abs(Z) <= Z_m then
      return math.max(NSrcFloor, (NSrc / math.sqrt(2. * math.pi * (sigSrc ^ 2))) *
         math.exp(-((lineLength - lineLengthSrc) ^ 2) / (2. * (sigSrc ^ 2))))
   else
      return 1.e-16
   end
end
srcTempIon           = function(t, xn)
   local lineLength = xn[3]
   local sigSrc = sigSrcIon
   local TSrc0, Tfloor = TSrc0Ion, TSrcFloorIon
   if math.abs(lineLength) <= 2. * sigSrc then
      return TSrc0
   else
      return Tfloor
   end
end

local dlineLengthEff = (lineLengthMax - lineLengthMin) / (numCellLineLength * (polyOrder + 1))
local kParMax        = math.pi / dlineLengthEff
local omega_H        = kParMax * vte / kperpRhos

log("\n")
log(string.format("  omega_H = kpar*vte/(kperp*rhos) = %e rad/s\n", omega_H))
log("\n")

local muMax_i, muMax_e = mi * ((3. * vti) ^ 2) / (2 * B_p), me * ((3. * vte) ^ 2) / (2 * B_p)
local tFinal           = 1e-20
local numFrames        = 1
-- Simulation App.
plasmaApp              = Plasma.App {
   logToFile          = true,
   --maxWallTime       = wallTime,                                        -- Wallclock time allocated for this job.
   tEnd               = tFinal,                                         -- End time.
   nFrame             = numFrames,                                      -- Number of output frames.
   lower              = { psiMin, -math.pi, lineLengthMin },            -- Configuration space lower left.
   upper              = { psiMax, math.pi, lineLengthMax },             -- Configuration space upper right.
   periodicDirs       = { 2 },                                          -- Periodic in theta.
   cells              = { numCellPsi, numCellTheta, numCellLineLength }, -- Configuration space cells.
   basis              = "serendipity",                                  -- One of "serendipity" or "maximal-order".
   polyOrder          = polyOrder,                                      -- Polynomial order.
   decompCuts         = { 1, 1, 1 },                                    -- Based on the number of cores available. Was {1,1,288}
   mapc2p             = function(xc)
      local psi, theta, lineLength = xc[1], xc[2], xc[3]                -- Field-aligned coordinates.
      local Z = Z_psiLineLength(psi, lineLength)                        -- Cylindrical axial coordinate Z.
      local R = R_psiZ(psi, Z)                                          -- Cylindrical radial coordinate R.
      -- Cartesian coordinates on plane perpendicular to Z axis.
      local X = R * math.cos(theta)
      local Y = R * math.sin(theta)
      return X, Y, Z
   end,

   timeStepper        = "rk3", -- One of "rk2" or "rk3".
   cflFrac            = 0.9,
   restartFrameEvery  = .05,
   calcIntQuantEvery  = 1. / (10. * numFrames), -- Aim at 10x more frequently than frames.
   parallelizeSpecies = true,
   -- Gyrokinetic electrons.
   elc                = Plasma.Species {
      charge = qe, mass = me,
      lower = { -3.75 * vte, 0 },
      upper = { 3.75 * vte, muMax_e },
      cells = { 96, 192 },
      -- Initial conditions.
      -- Will change when loaded from file.
      init = Plasma.MaxwellianProjection {
         --         fromFile = "elcIC.bp",
         density = function(t, xn)
            local psi, lineLength = xn[1], xn[3]
            local Z               = Z_psiLineLength(psi, lineLength) -- Cylindrical axial coordinate.
            local R               = R_psiZ(psi, Z)                   -- Cylindrical radial coordinate.
            local _, _, Bmag      = Bfield_psiZ(psi, Z)
            -- Density profile based on initial conditions
            if math.abs(Z) <= Z_bt then
               return n0 * math.sqrt((1. - ((R - R_bt) / alim) ^ 2) ^ alphaIC[1])
            elseif math.abs(Z) <= Z_m then
               return n0 * math.sqrt((1. - ((R - R_bt) / alim) ^ 2) ^ alphaIC[2])
            else
               return n_m * math.sqrt(Bmag / B_m)
            end
         end,
         driftSpeed = function(t, xn)
            local lineLength = xn[3]
            return cs_m * (lineLength / lineLength_m)
         end,
         temperature = function(t, xn)
            local psi, lineLength = xn[1], xn[3]
            local Z               = Z_psiLineLength(psi, lineLength) -- Cylindrical axial coordinate.
            local R               = R_psiZ(psi, Z)                   -- Cylindrical radial coordinate.
            local _, _, Bmag      = Bfield_psiZ(psi, Z)
            if math.abs(Z) <= Z_bt then
               return Te0 * math.sqrt((1. - ((R - R_bt) / alim) ^ 2) ^ alphaIC[1])
            elseif math.abs(Z) <= Z_m then
               return Te0 * math.sqrt((1. - ((R - R_bt) / alim) ^ 2) ^ alphaIC[2])
            else
               return Te_m * math.sqrt(Bmag / B_m)
            end
         end,
         --         density = "ion_gridDiagnostics_400.bp",
         --         driftSpeed = "ion_gridDiagnostics_400.bp",
         --         temperature = "ion_gridDiagnostics_400.bp",
         scaleWithSourcePower = false,
      },
      coll = Plasma.LBOCollisions {
         collideWith = { 'elc' },
         frequencies = { nuElc },
      },
      source = Plasma.Source {
         --         fromFile    = "elc_fSourceIC.bp",
         density     = srcDenElc,
         temperature = srcTempElc,
         diagnostics = { "twoFiles", "M0", "intM0", "intM2" },
      },
      --      polarizationDensityFactor = "ion_polarizationDensityFactor_0.bp",
      evolve = true,
      randomseed = randomseed,
      diagnostics = { "twoFiles", "M0", "Upar", "Temp", "Tperp", "Tpar", "intM0", "intM1", "intM2", "intKE", "intEnergy" },
      bcx = { Plasma.AbsorbBC {}, Plasma.ZeroFluxBC {} },
      bcz = {
         Plasma.SheathBC { diagnostics = { "twoFiles", "M0", "M1", "M2", "Upar", "Temp", "Energy", "intM0", "intM1",
            "intKE", "intEnergy" } },
         Plasma.SheathBC { diagnostics = { "twoFiles", "M0", "M1", "M2", "Upar", "Temp", "Energy", "intM0", "intM1",
            "intKE", "intEnergy" } } },
   },
   -- Gyrokinetic ions.
   ion                = Plasma.Species {
      charge = qi, mass = mi,
      lower = { -3.75 * vti, 0 },
      upper = { 3.75 * vti, muMax_i },
      cells = { 96, 192 },
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         --         fromFile = "ionIC.bp",
         density = function(t, xn)
            local psi, lineLength = xn[1], xn[3]
            local Z               = Z_psiLineLength(psi, lineLength) -- Cylindrical axial coordinate.
            local R               = R_psiZ(psi, Z)                   -- Cylindrical radial coordinate.
            local _, _, Bmag      = Bfield_psiZ(psi, Z)
            -- Density profile based on initial conditions
            if math.abs(Z) <= Z_bt then
               return n0 * math.sqrt((1. - ((R - R_bt) / alim) ^ 2) ^ alphaIC[1])
            elseif math.abs(Z) <= Z_m then
               return n0 * math.sqrt((1. - ((R - R_bt) / alim) ^ 2) ^ alphaIC[2])
            else
               return n_m * math.sqrt(Bmag / B_m)
            end
         end,
         driftSpeed = function(t, xn)
            local lineLength = xn[3]
            return cs_m * (lineLength / lineLength_m)
         end,
         temperature = function(t, xn)
            local psi, lineLength = xn[1], xn[3]
            local Z               = Z_psiLineLength(psi, lineLength) -- Cylindrical axial coordinate.
            local R               = R_psiZ(psi, Z)                   -- Cylindrical radial coordinate.
            local _, _, Bmag      = Bfield_psiZ(psi, Z)
            if math.abs(Z) <= Z_bt then
               return Ti0 * math.sqrt((1. - ((R - R_bt) / alim) ^ 2) ^ alphaIC[1])
            elseif math.abs(Z) <= Z_m then
               return Ti0 * math.sqrt((1. - ((R - R_bt) / alim) ^ 2) ^ alphaIC[2])
            else
               return Ti_m * math.sqrt(Bmag / B_m)
            end
         end,
         --         density = "ion_gridDiagnostics_400.bp",
         --         driftSpeed = "ion_gridDiagnostics_400.bp",
         --         temperature = "ion_gridDiagnostics_400.bp",
         scaleWithSourcePower = false,
      },
      coll = Plasma.LBOCollisions {
         collideWith = { 'ion' },
         frequencies = { nuIon },
      },
      source = Plasma.Source {
         --         fromFile    = "ion_fSourceIC.bp",
         density     = srcDenIon,
         temperature = srcTempIon,
         diagnostics = { "twoFiles", "M0", "intM0", "intM2" },
      },
      --      polarizationDensityFactor = "ion_polarizationDensityFactor_0.bp",
      evolve = true,
      randomseed = randomseed,
      diagnostics = { "twoFiles", "M0", "Upar", "Temp", "Tperp", "Tpar", "intM0", "intM1", "intM2", "intKE", "intEnergy" },
      bcx = { Plasma.AbsorbBC {}, Plasma.ZeroFluxBC {} },
      bcz = {
         Plasma.SheathBC { diagnostics = { "twoFiles", "M0", "M1", "M2", "Upar", "Temp", "Energy", "intM0", "intM1",
            "intKE", "intEnergy" } },
         Plasma.SheathBC { diagnostics = { "twoFiles", "M0", "M1", "M2", "Upar", "Temp", "Energy", "intM0", "intM1",
            "intKE", "intEnergy" } } },
   },

   -- Field solver.
   field              = Plasma.Field {
      -- Dirichlet in psi, Neuman in lineLength
      bcLowerPhi             = { { T = "D", V = 0.0 }, { T = "P", V = 0.0 } },
      bcUpperPhi             = { { T = "D", V = 0.0 }, { T = "P", V = 0.0 } },
      evolve                 = true, -- Evolve fields?
      uniformPolarization    = false,
      linearizedPolarization = false,
   },

   -- Magnetic geometry.
   funcField          = Plasma.Geometry {
      -- Background magnetic field.
      bmag   = function(t, xn)
         local psi, lineLength = xn[1], xn[3]
         local Z               = Z_psiLineLength(psi, lineLength)
         local _, _, Bmag      = Bfield_psiZ(psi, Z)
         return Bmag
      end,
      evolve = false, -- Geometry is not time-dependent.
      --      fromFile = "allGeoIC.bp",
   },

}
-- Run application.
plasmaApp:run()