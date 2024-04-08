-- Gkyl ------------------------------------------------------------------------
local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

gasGamma = 1.4 -- gas adiabatic constant

-- create app
eulerApp = Moments.App.new {
   logToFile = true,

   tEnd = 0.8, -- end time
   nFrame = 4, -- number of output frame
   lower = {0.0, 0.0}, -- lower left corner
   upper = {1.0, 1.0}, -- upper right corner
   cells = {200, 200}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   
   -- decomposition stuff
   decompCuts = {1, 1}, -- cuts in each direction

   -- electrons
   fluid = Moments.Species.new {
      charge = 0.0, mass = 1.0,
      equation = Euler.new { gasGamma = gasGamma },

      -- initial conditions
      init = function (t, xn)
	 -- Case 3 from Liska and Wendroff 2003 (see Table 4.3) these tables
	 -- store pressure, rho, u and v
	 local upLeft = {0.3, 0.5323, 1.206, 0.0}
	 local upRight = {1.5, 1.5, 0.0, 0.0}
	 local loLeft = {0.029, 0.138, 1.206, 1.206}
	 local loRight = {0.3, 0.5323, 0.0, 1.206}
	 local rho, u, v, pr      
	 local x, y = xn[1], xn[2]
	 local loc = 0.8
	 if y>loc then
	    if x<loc then
	       pr, rho, u, v = upLeft[1], upLeft[2], upLeft[3], upLeft[4]
	    else
	       pr, rho, u, v = upRight[1], upRight[2], upRight[3], upRight[4]
	    end
	 else
	    if x<loc then
	       pr, rho, u, v = loLeft[1], loLeft[2], loLeft[3], loLeft[4]
	    else
	       pr, rho, u, v = loRight[1], loRight[2], loRight[3], loRight[4]
	    end
	 end
	 return rho, rho*u, rho*v, 0.0, pr/(gasGamma-1) + 0.5*rho*(u^2+v^2)
      end,
      
      evolve = true, -- evolve species?

      bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
      bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
   },   
}
-- run application
eulerApp:run()
