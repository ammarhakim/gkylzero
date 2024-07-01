import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import postgkyl as pg

#simName='step_outboard_shaped_plate'
#simName='step_outboard_horizontal_plate'
#simName='step_outboard_fixed_z'
#simName='stepbry'
#simName='stepcore'

plt.figure()
simNames = ['wham']
for simName in simNames:
    dir = "./"
    data = pg.GData(dir+"xyz"+simName+"_nodes.gkyl")
    vals = data.getValues()
    X = vals[:,0,:,0]
    Y = vals[:,0,:,1]
    Z = vals[:,0,:,2]
    R=np.sqrt(X**2+Y**2)
    
    psid = pg.GData("wham_psi.gkyl")
    interp = pg.GInterpModal(psid,1,"ms")
    grid, psi = interp.interpolate()
    #pg.output.plot(psid, contour=True)
    
    for d in range(len(grid)):
        grid[d] = 0.5*(grid[d][:-1] + grid[d][1:])
    
    psisep= 0.0026
    psi_min = 0.001
    psi_max =  0.003
    npsi = 3
    dpsi = (psi_max-psi_min)/npsi
    
    clevels = np.linspace(psi_min, psi_max, npsi)
    plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="k")
    plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=np.r_[psisep], colors="r")
    
    #clevels = np.linspace(psi_min-3*dpsi, psi_max, npsi+4)
    #plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="r")
    
    plt.plot(R,Z,marker=".", color="k", linestyle="none")
    plt.scatter(R,Z, marker=".")
    segs1 = np.stack((R,Z), axis=2)
    segs2 = segs1.transpose(1,0,2)
    plt.gca().add_collection(LineCollection(segs1))
    plt.gca().add_collection(LineCollection(segs2))
    plt.grid()
    #plt.axis("tight")
    #plt.axis("image")
plt.savefig("step-dn.png", dpi=300)
plt.show()
