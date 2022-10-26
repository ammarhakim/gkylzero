-- Gkyl ------------------------------------------------------------------------
local Moments = require "Moments"

gasGamma = 2.0

-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 0.1, -- end time
   nFrame = 1, -- number of output frame
   lower = { 0.0 }, -- lower left corner
   upper = { 1.0 }, -- upper right corner
   cells = { 400 }, -- number of cells

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Moments.MHD {
	 gasGamma = gasGamma,
	 rp_type = "roe", -- one of "roe", "hlld", "lax"
	 divergenceConstraint = "none",
      },

      -- initial conditions
      init = function (t, xn)
	 local bx = 0.75
	 local rhol, rhor = 1.0, 0.125
	 local byl, byr = 1.0, -1.0
	 local pl, pr = 1.0, 0.1

	 local rho, by, p = rhor, byr, pr
	 if xn[1] < 0.5 then
	    rho, by, p = rhol, byl, pl
	 end
	 
	 return rho, 0.0, 0.0, 0.0, p/(gasGamma-1)+0.5*(bx^2+by^2), bx, by, 0.0
      end,
      
      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },
}
-- run application
eulerApp:run()
