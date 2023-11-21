import os
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    numRuns = 1

    for r in range(0, numRuns):
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
            "Internal Energy",
            "Gamma (C)",
            "Gamma (E)"
        ]

        # Loop through each time step
        for (j, file) in enumerate(datFiles):
            # Make subplots
            fig, axs = plt.subplots(3,2)

            # Load current timestep
            A = np.loadtxt(file)

            # Get the x-coords
            x = A[:,0]

            # Loop through each subplot/var
            for i, (ax, name) in enumerate(
                zip(
                axs.flatten(),
                primitiveVarNames
                )
                ):
                # Get curent var values
                if (name == "Internal Energy"):
                    var = A[:,i + 1]/A[:,1]
                else:
                    var = A[:,i + 1]

                # Plot the current var values
                ax.plot(x, var, "r--o",markerfacecolor="none")
                ax.set_xlabel("X", size=16)
                ax.set_ylabel("U", size=16)
                ax.set_title(name, size=18)
                ax.grid()
            
            # Show and then automatically close plot
            # Fake animation (improve later)
            if (j != len(datFiles) - 1):
                # plt.tight_layout()
                # plt.show(block=False)
                # plt.pause(0.01)
                # plt.close()
                pass
            else:
                plt.tight_layout()
                plt.savefig(f"../../Python/plots/plot{r + 1}.png", dpi=250)
                plt.show(block=True)
        os.chdir("../../Python")
