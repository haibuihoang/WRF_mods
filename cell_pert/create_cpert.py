from netCDF4 import Dataset
run_dir="./"
from wrf import getvar
from math import ceil
import numpy as np
import matplotlib.pyplot as plt

#This should be read from namelist
cell_pert=[0,0,0,1]
ndom = 4

vert_cpert = False  # If vert_cell, the the perturbation will be apply for every 8 grid cell in the vertical
#Fix option
zmin=400
zmax=1000
cell_sz=8

# The pertubation

def main():
    for i,cpert in enumerate(cell_pert):    
        print(f"Checking domain D{i+1:02}")
        if (cpert):
            wrfin = f"{run_dir}/wrfinput_d{i+1:02}"
            print("Creating cell_perturbation for ",wrfin)
            apply_cell_pert(wrfin)
            
            #Check
            #dset = Dataset(wrfin, 'r')
            #cell_pert_TH = dset['CELL_PERT_TH'][0,:,:,:]
            #plt.pcolormesh(cell_pert_TH[:,0,:],cmap="RdBu")

def apply_cell_pert(ifile):
    dset = Dataset(ifile, 'r+')
    cell_pert_TH = dset['CELL_PERT_TH'][0,:,:,:]
    h = getvar(dset,"z")[:,0,0].data
    kmax = np.argmin(np.abs(zmax-h))  # Only calculate form k=0 to kmax
    
    print(cell_pert_TH.shape,kmax)
    if (vert_cpert):
        N_cz = ceil(kmax/cell_sz)
        print(kmax,N_cz)
        for k in range(N_cz):
            k1,k2 = k*cell_sz, (k+1)*cell_sz
            temp_pert = gen_cell_turb_bdy(cell_pert_TH[k1,:,:]) 
            cell_pert_TH[k1:k2,:,:] = np.broadcast_to(temp_pert,cell_pert_TH[k1:k2,:,:].shape)
    else:
        for k in range(kmax+1):
            cell_pert_TH[k,:,:]=gen_cell_turb_bdy(cell_pert_TH[k,:,:]) 
        
    wgt = get_weight(h,zmin,zmax)
    for k in range(len(h)):
        cell_pert_TH[k,:,:] = cell_pert_TH[k,:,:]*wgt[k]
        
        
    dset['CELL_PERT_TH'][0,:,:,:] =  cell_pert_TH
    dset.close()
        
    
##Functions
def gen_cell_turb(perturb,pert_mag):
    pert_max=pert_mag
    pert_min=-pert_max
    dims=perturb.shape
    Nx,Ny = dims[1],dims[0]
    N_cx=ceil(Nx/cell_sz)
    N_cy=ceil(Ny/cell_sz)
    
    cell_perturb=pert_min + (pert_max-pert_min)*np.random.random((N_cy,N_cx)) 
    for i in range(N_cx):
        for j in range(N_cy):
            i1=i*cell_sz; i2=i1+cell_sz
            j1=j*cell_sz; j2=j1+cell_sz
            if (i2>Nx):
                i2=Nx
            if (j2>Ny):
                j2=Ny
            perturb[j1:j2,i1:i2]=cell_perturb[j,i]
    return perturb
    
def gen_cell_turb_bdy(perturb,pert_mag=0.5,bdy=24):
    dims=perturb.shape
    Nx,Ny = dims[1],dims[0]
    i1=j1=bdy
    j2=Ny-bdy
    i2=Nx-bdy
    print(perturb[:j1,:i2].dtype)
    perturb[:j1,:i2] = gen_cell_turb(perturb[:j1,:i2],pert_mag)
    perturb[:j2,i2:] = gen_cell_turb(perturb[:j2,i2:],pert_mag)
    perturb[j2:,i1:] = gen_cell_turb(perturb[j2:,i1:],pert_mag)
    perturb[j1:,:i1] = gen_cell_turb(perturb[j1:,:i1],pert_mag)
    return perturb

def get_weight(z,zmin,zmax):
    from math import pi
    kmin=np.argmin(np.abs(z-zmin))
    kmax=np.argmin(np.abs(z-zmax))
    wgt = np.cos(pi*(z-zmin)/(zmax-zmin)/2)**2
    wgt[0:kmin]=1
    wgt[kmax:]=0
    return wgt

# Run the main
main()
