import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

if __name__ == "__main__":
    numRuns = 22

    #for r in range(0, numRuns):
    # for r in range(0, numRuns):
    #for r in range(32, 33):
    #for r in range(25, 26):
    for r in range(41, 46):
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

        A = np.loadtxt(datFiles[0])
        x = A[:, 0]

        fig, axs = plt.subplots(4,1, figsize=(16,8))
        for i, (symbol, ax) in enumerate(zip(primitiveVarSymbols, axs.flatten())):
            if i != 3:
                ax.grid()
                ax.set_xlabel("X", size=16)
                ax.set_ylabel(symbol, size=16)
                ax.set_ylim([np.min(A[:, i + 1]) - 1.0, np.max(A[:, i + 1]) + 1.0])
            #ax.set_ylim([np.min(A[:, i + 1] - 5.0), np.max(A[:, i + 1]) + 25.0])
        lines = [axs[i].plot(x, A[:, i + 1], linestyle="--", markerfacecolor="none")[0] for i in range(0, 3)]
        div = np.zeros(len(A[:,2]))
        div[1:-1] = (A[2:, 2] - A[0:-2, 2])/(2 * 1/512)
        axs[3].set_ylim([np.min(div) - 1.0, np.max(div) + 1.0])
        lines.append(axs[3].plot(x,div)[0])
        points0 = [axs[0].plot(x_val, var_val, marker='o', markeredgecolor="blue", markerfacecolor="none")[0] for x_val, var_val in zip(x, A[:, 1])]
        points1 = [axs[1].plot(x_val, var_val, marker='o', markeredgecolor="blue", markerfacecolor="none")[0] for x_val, var_val in zip(x, A[:, 2])]
        points2 = [axs[2].plot(x_val, var_val, marker='o', markeredgecolor="blue", markerfacecolor="none")[0] for x_val, var_val in zip(x, A[:, 3])]
        pointss = [points0, points1, points2]


        plt.tight_layout()

        def update(frame):
            curr_file = datFiles[frame]
            A = np.loadtxt(curr_file)
            x = A[:, 0]
            preds = A[:,-1]

            div = np.zeros(len(A[:,2]))
            #print( (A[2:, 2] - A[0:-2, 2]))
            div[1:-1] = (A[2:, 2] - A[0:-2, 2])/(2 * 1/512)

            lines[3].set_ydata(div)

            for i, (line, points, ax) in enumerate(zip(lines, pointss, axs.flatten())):
                if i!=3:
                    var = A[:, (i + 1)]
                    line.set_ydata(var)
                    ax.set_ylim([np.min(var) - 1.0, np.max(var) + 1.0])
                    for point, var_val, pred in zip(points, A[:, (i+1)], preds):
                        if (pred == 1):
                            c = "red"
                        elif (pred == 2):
                            c = "pink"
                        else:
                            c = "blue"
                        point.set_ydata(var_val)
                        point.set_markeredgecolor(c)

            axs[3].set_ylim([np.min(div) - 5.0, 0])
            return lines[0], lines[1], lines[2], lines[3]



        # for i, (symbol, ax) in enumerate(zip(primitiveVarSymbols, axs.flatten())):
        #     var = A[:, (i + 1)]

        #     #ax.plot(x, var, "r--o",markerfacecolor="none")
        #     ax.plot(x, var, linestyle="--", markerfacecolor="none")
        #     for x_val, var_val, pred in zip(x, A[:, (i+1)], preds):
        #         if (pred == 1):
        #             c = "red"
        #         else:
        #             c = "blue"
        #         ax.plot(x_val, var_val, marker='o' ,markeredgecolor=c, markerfacecolor="none")

        #     ax.set_xlabel("X", size=16)
        #     ax.set_ylabel(symbol, size=16)
        #     #ax.set_title(name, size=18)
        #     ax.grid()


        ani = animation.FuncAnimation(fig=fig, func=update, frames=len(datFiles), interval=30)
        ani.save(f"../../Python/animations/anim{r+1}.mp4")
        os.chdir("../../Python")


