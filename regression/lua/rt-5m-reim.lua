-- Gkyl ------------------------------------------------------------------------
local Moments = require "Moments"

-- global parameters
charge = 10.0
gasGamma = 5./3.
elcCharge = -charge
ionCharge = charge
elcMass = 1/1836.2
ionMass = 1.0
lightSpeed = 1.0
epsilon0 = 1.0

local momentApp = Moments.App {

   tEnd = 10.0,
   nFrame = 1,
   lower = { 0.0 },
   upper = { 1.0 },
   cells = { 1000 },

   -- electrons
   elc = Moments.Species {
      charge = elcCharge, mass = elcMass,

      equation = Moments.Euler { gasGamma = gasGamma },
      
      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]

	 local sloc = 0.5 -- location of shock
	 if (x<sloc) then
	    return 1.0*elcMass/ionMass, 0.0, 0.0, 0.0, 5.0e-5/(gasGamma-1)
	 else
	    return 0.125*elcMass/ionMass, 0.0, 0.0, 0.0, 5.0e-6/(gasGamma-1)
	 end
      end,

      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },

   -- ions
   ion = Moments.Species {
      charge = ionCharge, mass = ionMass,

      equation = Moments.Euler { gasGamma = gasGamma },
      
      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]

	 local sloc = 0.5 -- location of shock
	 if (x<sloc) then
	    return 1.0, 0.0, 0.0, 0.0, 5.0e-5/(gasGamma-1)
	 else
	    return 0.125, 0.0, 0.0, 0.0, 5.0e-6/(gasGamma-1)
	 end
      end,

      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },

   field = Moments.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      
      init = function (t, xn)
	 local x, y = xn[1], xn[2]

	 local sloc = 0.5 -- location of shock
	 if (x<sloc) then
	    return 0.0, 0.0, 0.0, 0.75e-2, 0.0, 1.0e-2
	 else
	    return 0.0, 0.0, 0.0, 0.75e-2, 0.0, -1.0e-2
	 end
      end,
      
      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },

}
-- run application
momentApp:run()
