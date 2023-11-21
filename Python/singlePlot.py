import os
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    numRuns = 22

    #for r in range(0, numRuns):
    # for r in range(0, numRuns):
    for r in range(21, 22):
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

        # Set name of variables in order
        primitiveVarNames = [
            "Density",
            "Velocity",
            "Pressure",
            #"Internal Energy",
            #"Gamma (C)",
            #"Gamma (E)"
        ]

        primitiveVarSymbols = [
            r"$\rho$",
            r"U",
            r"P"
        ]

        finalFile = datFiles[-1]

        os.system(f"cp ./{finalFile} ../out/out{r + 1}.txt ")

        A = np.loadtxt(finalFile)
        x = A[:, 0]

        fig, axs = plt.subplots(3,1)

        for i, (symbol, ax) in enumerate(zip(primitiveVarSymbols, axs.flatten())):
            var = A[:, (i + 1)]

            #ax.plot(x, var, "r--o",markerfacecolor="none")
            ax.plot(x, var, markerfacecolor="none")
            ax.set_xlabel("X", size=16)
            ax.set_ylabel(symbol, size=16)
            #ax.set_title(name, size=18)
            ax.grid()

        plt.tight_layout()
        plt.savefig(f"../../Python/plots/plot{r + 1}.png", dpi=250)
        #plt.show()
        plt.close()
        os.chdir("../../Python")
