-- Gkyl ------------------------------------------------------------------------
local Moments = require "Moments"

gasGamma = 1.4

-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 0.1, -- end time
   nFrame = 1, -- number of output frame
   lower = { 0.0 }, -- lower left corner
   upper = { 1.0 }, -- upper right corner
   cells = { 100 }, -- number of cells

   scheme_type = "mp",
   mp_recon = "u5",

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Moments.Euler {
	 gasGamma = gasGamma,
	 rp_type = "lax", -- one of "roe", "hllc", "lax"
      },

      -- initial conditions
      init = function (t, xn)
	 local rhol, ul, pl = 3.0, 0.0, 3.0
	 local rhor, ur, pr = 1.0, 0.0, 1.0

	 local rho, u, p = rhor, ur, pr
	 if xn[1] < 0.5 then
	    rho, u, p = rhol, ul, pl
	 end
	 
	 return rho, rho*u, 0.0, 0.0, 0.5*rho*u^2 + p/(gasGamma-1)
      end,
      
      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },
}
-- run application
eulerApp:run()
