local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

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

OmegaPe0 = math.sqrt(n0 * elcCharge^2 / (epsilon0 * elcMass))
de0 = lightSpeed / OmegaPe0
OmegaPi0 = math.sqrt(n0 * ionCharge^2 / (epsilon0 * ionMass))
di0 = lightSpeed / OmegaPi0
lambda = lambdaOverDi0 * di0

-- domain size
Lx = 25.6 * di0
Ly = 12.8 * di0

momentApp = Moments.App.new {

   tEnd = 25.0/OmegaCi0,
   nFrame = 1,
   lower = {-Lx/2, -Ly/2},
   upper = {Lx/2, Ly/2},
   cells = {2*64, 2*32},
   cflFrac = 1.0,

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   elc = Moments.Species.new {
      charge = elcCharge, mass = elcMass,

      equation = Euler.new { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]

	 local TeFrac = 1.0/(1.0 + TiOverTe)
	 local sech2 = (1.0/math.cosh(y/lambda))^2
	 local n = n0*(sech2 + nbOverN0)
	 local Jz = -(B0/lambda)*sech2
	 local Ttotal = plasmaBeta*(B0*B0)/2.0/n0

	 local rhoe = n*elcMass
	 local ezmom = (elcMass/elcCharge)*Jz*TeFrac
	 local ere = n*Ttotal*TeFrac/(gasGamma-1) + 0.5*ezmom*ezmom/rhoe
	 
	 return rhoe, 0.0, 0.0, ezmom, ere
      end,
      evolve = true, -- evolve species?
      bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall },
   },

   -- ions
   ion = Moments.Species.new {
      charge = ionCharge, mass = ionMass,

      equation = Euler.new { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]

	 local TiFrac = TiOverTe/(1.0 + TiOverTe)
	 local sech2 = (1.0/math.cosh(y/lambda))^2
	 local n = n0*(sech2 + nbOverN0)
	 local Jz = -(B0/lambda)*sech2
	 local Ttotal = plasmaBeta*(B0*B0)/2.0/n0

	 local rhoi = n*ionMass
	 local izmom = (ionMass/ionCharge)*Jz*TiFrac
	 local eri = n*Ttotal*TiFrac/(gasGamma-1) + 0.5*izmom*izmom/rhoi
	 
	 return rhoi, 0.0, 0.0, izmom, eri
      end,
      evolve = true, -- evolve species?
      bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall },
   },

   field = Moments.Field.new {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local x, y = xn[1], xn[2]

	 local Bxb = B0 * math.tanh(y / lambda) 
	 local Bx = Bxb - psi0 *(math.pi / Ly) * math.cos(2 * math.pi * x / Lx) * math.sin(math.pi * y / Ly) 
	 local By = psi0 * (2 * math.pi / Lx) * math.sin(2 * math.pi * x / Lx) * math.cos(math.pi * y / Ly)
	 local Bz = 0.0

	 return 0.0, 0.0, 0.0, Bx, By, Bz
      end,
      evolve = true, -- evolve field?
      bcy = { G0.FieldBc.bcWall, G0.FieldBc.bcWall },
   },

}

-- check each Lua-C interface method

field_ic_status = momentApp:apply_ic_field()
print("field_ic_status", field_ic_status)

elc_ic_status = momentApp:apply_ic_species("elc")
print("elc_ic_status", elc_ic_status)

ion_ic_status = momentApp:apply_ic_species("ion")
print("ion_ic_status", ion_ic_status)

field_io_status = momentApp:write_field(0.0, 0)
print("field_io_status", field_io_status)

elc_io_status = momentApp:write_species("elc", 0.0, 0)
print("elc_io_status", elc_io_status)

ion_io_status = momentApp:write_species("ion", 0.0, 0)
print("ion_io_status", ion_io_status)

ic_status = momentApp:apply_ic()
print("ic_status", ic_status)

io_status = momentApp:write(1.0, 1)
print("io_status", io_status)

max_dt = momentApp:max_dt()
print("max_dt", max_dt)

field_energy_status = momentApp:calc_field_energy(0.0)
print("field_energy_status", field_energy_status)

int_mom_status = momentApp:calc_integrated_mom(0.0)
print("int_mom_status", int_mom_status)

write_field_energy_status = momentApp:write_field_energy()
print("write_field_energy_status", write_field_energy_status)

write_int_mom_status = momentApp:write_integrated_mom()
print("write_int_mom_status", write_int_mom_status)

stat_write_status = momentApp:stat_write()
print("stat_write_status", stat_write_status)

update_status = momentApp:update(0.1)
success = update_status and 1 or 0
io.write( string.format("update_status = { success = %d, dt_actual = %g, dt_suggested = %g }\n", 
  success, update_status.dt_actual, update_status.dt_suggested)
)

nghost = momentApp:nghost()
io.write( string.format("nghost = { %d, %d, %d }\n", nghost[1], nghost[2], nghost[3]) )

nfield_energy = momentApp:field_energy_ndiag()
print("Number of field-energy diagonostics are", nfield_energy)
