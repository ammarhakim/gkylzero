-- Gkyl ------------------------------------------------------------------------
local Moments = require "Moments"

gasGamma = 1.4 -- gas adiabatic constant

-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 0.1, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.25, 0.0}, -- lower left corner
   upper = {1.25, 2*math.pi}, -- upper right corner
   cells = {64, 64*6}, -- number of cells
   cflFrac = 0.9,

   periodicDirs = { 2 },

   mapc2p = function(t, xn)
      local r, th = xn[1], xn[2]
      return r*math.cos(th), r*math.sin(th)
   end,

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Moments.Euler { gasGamma = gasGamma },

      init = function (t, xn)
	 local r = xn[1]
	 local rhol, ul, pl = 3.0, 0.0, 3.0
	 local rhor, ur, pr = 1.0, 0.0, 1.0
	 local sloc = 0.5*(0.25+1.25)

	 local rho, u, pr = rhor, ur, pr
	 if r<sloc then
	    rho, u, pr = rhol, ul, pl
	 end

	 return rho, rho*u, 0.0, 0.0, pr/(gasGamma-1) + 0.5*rho*u^2
      end,

      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
      
   },
}
-- run application
eulerApp:run()
