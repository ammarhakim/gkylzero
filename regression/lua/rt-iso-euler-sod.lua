-- Gkyl ------------------------------------------------------------------------
local Moments = require "Moments"

cs = 1.0 -- thermal speed

-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 0.1, -- end time
   nFrame = 1, -- number of output frame
   lower = { 0.0 }, -- lower left corner
   upper = { 1.0 }, -- upper right corner
   cells = { 100 }, -- number of cells

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Moments.IsoEuler { vthermal = cs },
      
      -- initial conditions
      init = function (t, xn)
	 local rhol, ul = 3.0, 0.5
	 local rhor, ur = 1.0, 0.0

	 local rho, u = rhor, ur
	 if xn[1] < 0.5 then
	    rho, u = rhol, ul
	 end
	 
	 return rho, rho*u, 0.0, 0.0
      end,
      
      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },
}
-- run application
eulerApp:run()
