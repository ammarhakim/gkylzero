import postgkyl as pg
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.integrate as integrate
from scipy.interpolate import pchip_interpolate
import scipy.optimize as sco
import fortranformat as ff
from datetime import date
from typing import Any, Generator, Iterable, List, TextIO, Union

B0 = 0.55
R0 = 2.5
B_axis = B0
def psi_f(R, Z):
    x = R/R0
    y = Z/R0
    return - (0.00373804283369698*x**6*np.log(x) - 0.00574955335438162*x**6 - 0.0448565140043638*x**4*y**2*np.log(x) + 0.0503044260840945*x**4*y**2 + 0.017623348727471*x**4*np.log(x) + 0.0956643504553684*x**4 + 0.0299043426695758*x**2*y**4*np.log(x) - 0.0160920841654771*x**2*y**4 - 0.070493394909884*x**2*y**2*np.log(x) + 0.0644725519961134*x**2*y**2 - 7.00898484785135e-5*x**2*np.log(x) - 0.303766642191745*x**2 - 0.00199362284463839*y**6 + 0.0117488991516473*y**4 + 7.00898484785135e-5*y**2 + 0.0145368720253973)


#RZ box
NW = 257
NH = 257
RMIN,RMAX = 0.1*R0, 2.4*R0
ZMIN,ZMAX = -2.4*R0, 2.4*R0
RDIM = RMAX - RMIN
ZDIM = ZMAX - ZMIN
RLEFT=RMIN
ZMID = (ZMAX+ZMIN)/2.0
RMAXIS = 0.0
ZMAXIS = 0.0
SIMAG = psi_f(R0,0) 
SIBRY = 0.0
NPSI = NW #don't write
RCENTR = R0
BCENTR = B0
CURRENT = 0


#Solve GS in RZ coords
Rgrid = np.linspace(RMIN,RMAX,NW)
Zgrid = np.linspace(ZMIN,ZMAX,NH)

#rthetagrid = np.zeros((len(Rgrid),len(Zgrid),2))
psiRZ = np.zeros((len(Rgrid),len(Zgrid)))
for i,Ri in enumerate(Rgrid):
    for j,Zj in enumerate(Zgrid):
        psiRZ[i,j] = psi_f(Ri,Zj)

plt.figure()
plt.contour(Rgrid,Zgrid, psiRZ.T, levels = np.linspace(-1,0.3,12))
plt.colorbar()
plt.title("Psi calculated in RZ coords")
plt.show()




#PSI quantities
PSIGRID = np.linspace(SIMAG, SIBRY,NPSI)
FPOL = (B0*R0/Rgrid)*Rgrid # F = RBphi
FFPRIM = np.repeat(0.0, NPSI)
PPRIME = np.repeat(-1e-6,NPSI)
PRES = integrate.cumtrapz(PPRIME,PSIGRID,initial=0)
PSIZR = psiRZ.T


QPSI = np.zeros(NPSI)
for i in range(NPSI):
    QPSI[i] = 0



writeList = [NW, NH,                                        #3i4
             RDIM, ZDIM, RCENTR, RLEFT, ZMID,               #5E16.9
             RMAXIS, ZMAXIS, SIMAG, SIBRY, BCENTR,          #5e16.9
             CURRENT, SIMAG, 0, RMAXIS, 0,                  #5E16.9
             ZMAXIS, 0, SIBRY, 0, 0,                        #5E16.9
             FPOL, PRES, FFPRIM, PPRIME,                    #5E16.9
             PSIZR, QPSI]                                   #5E16.9

#Header stuff
header_fmt = "(a48,3i4)"
label='AAAAAA'
creation_date = date.today().strftime("%d/%m/%Y")
shot = int(0)
time = int(0)
shot_str = f"# {shot:d}"
time_str = f"  {time:d}ms"
comment = f"{label:11}{creation_date:10s}   {shot_str:>8s}{time_str:16s}"
def write_line(data: Iterable[Any], fh: TextIO, fmt: str) -> None:
    r"""
    Writes to a Fortran formatted ASCII data file. The file handle will be left on a
    newline.

    Parameters
    ---------
    data:
        The data to write.
    fh:
        File handle. Should be in a text write mode, i.e. ``open(filename, "w")``.
    fmt:
        A Fortran IO format string, such as ``'(6a8,3i3)'``.
    """
    fh.write(ff.FortranRecordWriter(fmt).write(data))
    fh.write("\n")
#write_line((comment, 3, NW, NH), fh, header_fmt) #3 is idum


#Now write the EFIT FILE
with open('cerfon.geqdsk','w',newline='') as f:
    #Comment, NW,NH
    write_line((comment, 3, NW, NH), f, header_fmt) #3 is idum
    # rdim,zdim,rcentr,rleft,zmid
    for i in range(2,7):
        f.write('%16.9E'%writeList[i])
    f.write('\n')
    # rmaxis,zmaxis,simag,sibry,bcentr
    for i in range(7,12):
        f.write('%16.9E'%writeList[i])
    f.write('\n')
    # current, simag, xdum, rmaxis, xdum
    for i in range(12,17):
        f.write('%16.9E'%writeList[i])
    f.write('\n')
    #zmaxis,xdum,sibry,xdum,xdum
    for i in range(17,22):
         f.write('%16.9E'%writeList[i])
    f.write('\n')
    #FPOL,PRES,FFPRIM,PPRIME
    for i in range(22,26):
        count=0
        for j in range(0,NW):
            f.write('%16.9E'%writeList[i][j])
            count = count+1
            if count==5:
                f.write('\n')
                count=0
    #PSIZR
    for i in range(26,27):
        count = 0
        for j in range(0,NH):
            for k in range(0,NW):
                f.write('%16.9E'%writeList[i][j][k])
                count = count+1
                if count==5:
                    f.write('\n')
                    count=0
    #QPSI
    for i in range(27,28):
        count=0
        for j in range(0,NW):
            f.write('%16.9E'%writeList[i][j])
            count = count+1
            if count==5:
                f.write('\n')
                count=0




plt.show()
