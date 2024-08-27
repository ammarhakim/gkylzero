import requests
import os
import shutil
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm

adas_data_dir = "."

# Adapted from Francesco Sciortino Copyright (c) 2021
class adas_adf11:

    def __init__(self,filename):
    #def __init__(self):     
        self.filepath = adas_data_dir + os.sep + filename
        self.filename = filename
        self.file_type = self.filename[:3]

        try:
            self.imp = self.filename.split("_")[1].split(".")[0]
        except:
            self.imp = None
        
        #get data
        self.load()
        #convert to lists to import to C
        self.getClist()
        
    def getClist(self):
        self.clogNe = self.logNe.tolist()
        self.clogT = self.logT.tolist()
        self.clogdata = self.logdata.ravel().tolist()
                
    #Load data
    def load(self):
        loc = get_adas_file_loc(self.filename)
        with open(loc) as f:
            header = f.readline()
            self.n_ion, self.n_ne, self.n_T = numpy.int_(header.split()[:3])
            details = " ".join(header.split()[3:])

            f.readline() #skip line of ---
            line = f.readline()
            #metastable resolved file
            if all([a.isdigit() for a in line.split()]):
                self.metastables = numpy.int_(line.split())
                f.readline() # skip empty line
                line = f.readline()
            else:
                self.metastables = numpy.ones(self.n_ion + 1, dtype = int)

            logNe = []
            while len(logNe) <self.n_ne:
                logNe += [float(n) for n in line.split()]
                line = f.readline()

            logT = []
            while len(logT) < self.n_T:
                logT += [float(t) for t in line.split()]
                line = f.readline()

            subheader = line

            ldata, self.Z, self.MGRD, self.MPRT=[],[],[],[]
            ind = 0
            while True:
                ind += 1

                try:#if file includs metastables states
                    iprt, igrd, typ, z = subheader.split("/")[1:5]
                    self.Z.append(int(z.split("=")[1]))
                    self.MGRD.append(int(igrd.split("=")[1]))
                    self.MPRT.append(int(iprt.split("=")[1]))
                except:
                    self.Z.append(ind + 1)
                    self.MGRD.append(1)
                    self.MPRT.append(1)
                    
                drcofd = [] #log10(generalized radiative coefficients)
                while len(drcofd) < self.n_ne * self.n_T:
                    line = f.readline()
                    drcofd += [float(L) for L in line.split()]

                ldata.append(numpy.array(drcofd).reshape(self.n_T,self.n_ne))

                subheader = f.readline().replace("-", " ")
                #end of the file
                if len(subheader) == 0 or subheader.isspace() or subheader[0] == "C":
                    break
                             
        self.logNe = numpy.array(logNe)
        self.logT = numpy.array(logT)
        self.logdata = numpy.array(ldata)

        self.meta_ind = list(zip(self.Z, self.MGRD, self.MPRT))

        def data_for_charge_state(self, charge_state):
            return self.logdata[int(charge_state),:,:]


    def plot(self, fig=None, axes=None):
        """Plot data from input ADAS file. If provided, the arguments allow users to overplot
        and compare data from multiple files.

        Parameters
        ----------
        fig : matplotlib Figure object
            If provided, add specification as to which ADAS file is being plotted.
        axes : matplotlib Axes object (or equivalent)
            If provided, plot on these axes. Note that this typically needs to be a set of axes
            for each plotted charge state. Users may want to call this function once first to get
            some axes, and then pass those same axes to a second call for another file to compare with.
        """

        # settings for plotting
        self.ncol = numpy.ceil(numpy.sqrt(len(self.Z)))
        self.nrow = numpy.ceil(len(self.Z) / self.ncol)

        if fig is None or axes is None:
            fig, axes = plt.subplots(
                int(self.ncol), int(self.nrow), sharex=True, sharey=True
            )

        axes = numpy.atleast_2d(axes)
        colormap = cm.rainbow
        colors = cm.rainbow(numpy.linspace(0, 1, len(self.logNe)))

        if fig is not None:
            fig.suptitle(self.filename + "  " + self.file_type)

        for i, ax in enumerate(axes.flatten()):
            if i >= len(self.Z):
                break
            if all(self.logdata[i].std(1) == 0):  # independent of density
                ax.plot(self.logT, self.logdata[i, :, 0])
            else:
                ax.set_prop_cycle("color", colors)
                ax.plot(self.logT, self.logdata[i])
                ax.text(
                    0.1,
                    0.8,
                    "$n_e = 10^{%.0f-%.0f}\mathrm{[cm^{-3}]}$"
                    % (self.logNe[0], self.logNe[-1]),
                    horizontalalignment="left",
                    transform=ax.transAxes,
                )

            ax.grid(True)

            if self.file_type != "brs":
                charge = self.Z[i]
                meta = self.MPRT[i], self.MGRD[i]
                if self.file_type in ["scd", "prs", "ccd", "prb", "qcd"]:
                    charge -= 1
                title = self.imp + "$^{%d\!+}$" % charge
                if any(self.metastables > 1):
                    title += str(meta)
                ax.set_title(title)

        for ax in axes[-1]:
            ax.set_xlabel("$\log\ T_e\ \mathrm{[eV]}$")
        for ax in axes[:, 0]:
            if self.file_type in ["scd", "acd", "ccd"]:
                ax.set_ylabel("$\log(" + self.file_type + ")\ \mathrm{[cm^3/s]}$")
            elif self.file_type in ["prb", "plt", "prc", "pls", "brs", "prs"]:
                ax.set_ylabel("$\log(" + self.file_type + ")\ \mathrm{[W\cdot cm^3]}$")        
        plt.show()
        


#Return full path for specific file and download if necessary
def get_adas_file_loc(filename):
    #if adas_data_dir doesn't exist - create it
    
    if filename == "none":
        #Don't load a file
        return
    elif os.path.exists(adas_data_dir + os.sep + filename):
        #File is already on system
        return adas_data_dir + os.sep + filename
    elif os.path.exists(filename):
        #filename is full path
        return filename
    else:
        loc = adas_data_dir + os.sep + filename
        fetch_file(filename,loc)
        return loc

    
#Fetch file from open.adas.ac.uk    
def fetch_file(filename,loc):
    url = "https://open.adas.ac.uk/download/adf11/"
    
    filename_mod = filename.split("_")[0] + "/" + filename
    
    r = requests.get(url + "/" + filename_mod)
    
    if(len(r.text)) < 1000:
        raise ValueError(f'Could not fetch {filename} from ADAS!')
    
    with open(loc, "wb") as f:
        f.write(r.content)
    

# 3d array is (Z, T, Ne)

def write_files(elem_array, file_num):

    for i in range(len(elem_array)):
        name = elem_array[i][0]
        num = file_num[elem_array[i][1]]
        zmax = elem_array[i][2]
        #print(name, num)
        ioniz_np = []
        recomb_np = []
        for zi in range(0,zmax):
            ioniz = adas_adf11("scd%s_%s.dat"%(num,name))
            ioniz_dat = ioniz.logdata[zi,:,:]-6.0
            ioniz_flat = numpy.ndarray.flatten(ioniz_dat)
            ioniz_np.append(ioniz_flat)
            
            recomb = adas_adf11("acd%s_%s.dat"%(num,name))
            recomb_dat = recomb.logdata[zi,:,:]-6.0
            recomb_flat = numpy.ndarray.flatten(recomb_dat)
            recomb_np.append(recomb_flat)
 
        #print('2d shape', numpy.shape(recomb_dat))
        ioniz_np = numpy.array(ioniz_np)
        recomb_np = numpy.array(recomb_np)
        ioniz_np.tofile("ioniz_%s.npy"%name)
        recomb_np.tofile("recomb_%s.npy"%name)
        ioniz.logT.tofile("logT_%s.npy"%name)
        ioniz.logNe.tofile("logN_%s.npy"%name)

file_num = ['12','96','89']
elem = [['h',0,1],['he',1,2],['li',1,3],['be',1,4],['b',2,5],['c',1,6],['n',1,7],['o',1,8],['al',2,13],['ar',2,18]]

write_files(elem, file_num)
