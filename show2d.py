import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":
    fname="tmax2.out";
    fname="amax2.out";
    fp=open(fname,"r");

    dat=fp.readline()

    dat=list(map(int,dat.strip().split(",")))

    Nx=dat[0];
    Ny=dat[1];

    amp=fp.read();
    amp=amp.strip();
    amp=np.array(amp.split("\n"));
    amp=amp.astype(np.float)
    amp=np.reshape(amp,[Nx,Ny])
    amp=amp.transpose();

    fig,ax=plt.subplots(1,1)
    ext=[-10,10,0,30]
    #ax.imshow(amp,aspect="equal",cmap="jet",vmin=-0.3,vmax=0.3,origin="lower",extent=ext,interpolation="bilinear")
    ax.imshow(amp,aspect="equal",cmap="jet",vmin=0.70,vmax=1.0,origin="lower",extent=ext)
    plt.show()
    

