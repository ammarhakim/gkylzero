-- Gkyl ------------------------------------------------------------------------
local Moments = require "Moments"

gasGamma = 5.0/3.0

-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 4.0, -- end time
   nFrame = 4, -- number of output frame
   lower = { 0.0 }, -- lower left corner
   upper = { 2.0 }, -- upper right corner
   cells = { 32 }, -- number of cells
   periodicDirs = { 1 },

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Moments.Euler { gasGamma = gasGamma },
      limiter = "no-limiter",
      
      -- initial conditions
      init = function (t, xn)
	 local x = xn[1]
	 local rho = 1.0 + 0.2*math.sin(math.pi*x)
	 local u, pr = 1.0, 1.0
	 return rho, rho*u, 0.0, 0.0, pr/(gasGamma-1) + 0.5*rho*u^2
      end,
   },
}
-- run application
eulerApp:run()
