import os
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    As = []
    for r in range(25, 29):
        # Clean the current output
        os.chdir("../Fortran/data/")
        os.system("rm *.dat")

        # Get the current input
        os.chdir("../inputs")
        currFile = f"input{r + 1}.txt"
        os.system("rm ../slug.init")
        os.system(f"cp ./{currFile} ../slug.init")

        # Got to the data folder
        os.chdir("../")
        os.system("./slugEuler1d")
        os.chdir("./data/")
        
        # Get a list of all files in Fortran dir
        files = [thing for thing in os.listdir() if os.path.isfile(thing)]

        # Get the dat files
        datFiles = np.array([file for file in files if file[-3:] == "dat"])


        # Get dat file numbers
        N = len(datFiles)
        nums = np.zeros(N)

        # Get step number from each file name
        for i, fileName in enumerate(datFiles):
            parts = fileName.split("_")
            keyPart = parts[3][:-4]
            num = int(keyPart)
            nums[i] = num
        
        # Sort the file names according to their number
        order = np.argsort(nums)
        datFiles = datFiles[order]

        finalFile = datFiles[-1]

        os.system(f"cp ./{finalFile} ../out/out{r + 1}.txt ")

        A = np.loadtxt(finalFile)
        As.append(A)

        plt.tight_layout()
        plt.savefig(f"../../Python/plots/plot{r + 1}.png", dpi=250)
        #plt.show()
        plt.close()
        os.chdir("../../Python")

    labels = [
        "NN (FOG + PLM)",
        "FOG",
        "PLM",
        "PPM",
    ]

    colors = [
        "orange",
        "green",
        "purple",
        "black"
    ]

    fig, ax = plt.subplots(figsize = (12,8))
    A = As[0]
    x = A[:, 0]
    dens = A[:, 1]
    #dens = A[:, 2]
    preds = A[:, -1]

    ax.plot(x, dens, linestyle = "--", label = labels[0])
    for xv, d, p in zip(x, dens, preds):
        if (p == 1):
            c = "red"
        elif (p == 2):
            c = "pink"
        else:
            c = "blue"
        ax.plot(xv, d, markerfacecolor = "none", markeredgecolor=c, marker='o')
    
    for i in range(1, 4):
        A = As[i]
        x = A[:, 0]
        dens = A[:, 1]
        #dens = A[:, 2]
        ax.plot(x, dens, linestyle="--", marker='o', markerfacecolor="none", label=labels[i], color = colors[i - 1])
    
    ax.grid()
    ax.set_xlabel("x",size=14)
    ax.set_ylabel(r"$\rho$",size=14)
    ax.legend()

    plt.savefig("./plots/blastCompareComb.png",dpi=200)

    fig, axs = plt.subplots(2, 2, figsize = (12,8))
    A = As[0]
    x = A[:, 0]
    dens = A[:, 1]
    #dens = A[:, 2]
    preds = A[:, -1]

    ax = axs.flatten()[0]
    ax.plot(x, dens, linestyle = "--", label = labels[0])
    for xv, d, p in zip(x, dens, preds):
        if (p == 1):
            c = "red"
        elif (p == 2):
            c = "pink"
        else:
            c = "blue"
        ax.plot(xv, d, markerfacecolor = "none", markeredgecolor=c, marker='o')
    ax.grid()
    ax.legend()
    ax.set_xlabel("x",size=14)
    ax.set_ylabel(r"$\rho$",size=14)
    
    for i in range(1, 4):
        ax = axs.flatten()[i]
        A = As[i]
        x = A[:, 0]
        dens = A[:, 1]
        #dens = A[:, 2]
        ax.plot(x, dens, linestyle="--", marker='o', markerfacecolor="none", label=labels[i], color = colors[i - 1])
        ax.grid()
        ax.set_xlabel("x",size=14)
        ax.set_ylabel(r"$\rho$",size=14)
        ax.legend()
    

    plt.savefig("./plots/blastCompareSep.png",dpi=200)

