-- Gkyl ------------------------------------------------------------------------
local Moments = require "Moments"

-- physical parameters
gasGamma = 5./3.
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25.0
lightSpeed = 1.0
epsilon0 = 1.0
mu0 = 1.0

n0 = 1.0
VAe = 0.5
plasmaBeta = 1.0
lambdaOverDi0 = 0.5
TiOverTe = 5.0
nbOverN0 = 0.2
pert = 0.1
Valf = VAe*math.sqrt(elcMass/ionMass)
B0 = Valf*math.sqrt(n0*ionMass)
OmegaCi0 = ionCharge*B0/ionMass
psi0 = pert*B0

OmegaPe0 = math.sqrt(n0*elcCharge^2/(epsilon0*elcMass))
de0 = lightSpeed/OmegaPe0
OmegaPi0 = math.sqrt(n0*ionCharge^2/(epsilon0*ionMass))
di0 = lightSpeed/OmegaPi0
lambda = lambdaOverDi0*di0

-- domain size
R0 = 10*di0 -- inner-radius
R1 = 25*di0 -- outer-radius

deltaPhi = 0.1 -- potential difference phi(R1)-phi(R0)

function radialE(r)
   return -deltaPhi/math.log(R1/R0)/r
end

-- create app
momentApp = Moments.App {
   
   tEnd = 100.0/OmegaPi0, -- end time
   nFrame = 20, -- number of output frame
   lower = {R0, -2*math.pi/8}, -- lower left corner
   upper = {R1, 2*math.pi/8}, -- upper right corner
   cells = {64, 64}, -- number of cells
   cflFrac = 0.9,

   --periodicDirs = { 2 },

   mapc2p = function(t, xn)
      local r, th = xn[1], xn[2]
      return r*math.cos(th), r*math.sin(th)
   end,

   -- electrons
   elc = Moments.Species {
      charge = elcCharge, mass = elcMass,

      equation = Moments.Euler { gasGamma = gasGamma },
      
      -- initial conditions
      init = function (t, xn)
	 local r, theta = xn[1], xn[2]
	 local Er = radialE(r)	 
	 local TeFrac = 1.0/(1.0 + TiOverTe)
	 local n = n0
	 local ux = Er/B0*math.sin(theta)
	 local uy = -Er/B0*math.cos(theta)	 
	 local Ttotal = plasmaBeta*(B0*B0)/2.0/n0

	 local rhoe = n*elcMass
	 local ere = 0.5*rhoe*(ux^2 + uy^2) + n*Ttotal*TeFrac/(gasGamma-1)
	 
	 return rhoe, rhoe*ux, rhoe*uy, 0.0, ere
      end,

      bcx = { Moments.Species.bcWall, Moments.Species.bcWall },
      bcy = { Moments.Species.bcWedge, Moments.Species.bcWedge },
   },

   -- ions
   ion = Moments.Species {
      charge = ionCharge, mass = ionMass,

      equation = Moments.Euler { gasGamma = gasGamma },
      
      -- initial conditions
      init = function (t, xn)
	 local r, theta = xn[1], xn[2]
	 local Er = radialE(r)
	 local TiFrac = TiOverTe/(1.0 + TiOverTe)
	 local n = n0
	 local ux = Er/B0*math.sin(theta)
	 local uy = -Er/B0*math.cos(theta)
	 local Ttotal = plasmaBeta*(B0*B0)/2.0/n0

	 local rhoi = n*ionMass
	 local eri = 0.5*rhoi*(ux^2 + uy^2) + n*Ttotal*TiFrac/(gasGamma-1)
	 
	 return rhoi, rhoi*ux, rhoi*uy, 0.0, eri
      end,
      
      bcx = { Moments.Species.bcWall, Moments.Species.bcWall },
      bcy = { Moments.Species.bcWedge, Moments.Species.bcWedge },      
   },

   field = Moments.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      
      init = function (t, xn)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,

      is_ext_em_static = true,
      ext_em_func = function (t, xn)
	 local r, theta = xn[1], xn[2]
	 local Er = radialE(r)
	 local Ex = Er*math.cos(theta)
	 local Ey = Er*math.sin(theta)
	 return Ex, Ey, 0.0, 0.0, 0.0, B0
      end,
      
      bcx = { Moments.Field.bcReflect, Moments.Field.bcReflect },
      bcy = { Moments.Field.bcWedge, Moments.Field.bcWedge },
   },
}
-- run application
momentApp:run()
