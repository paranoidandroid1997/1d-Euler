import os 
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    numFiles = 21
    data = []

    for i in range(1, numFiles + 1):
        currFile = f"../Fortran/out/out{i}.txt"
        data.append(np.loadtxt(currFile))

    # First Plot
    fig = plt.figure(constrained_layout=True, figsize=(16,8))
    fig.suptitle(
        r"Comparison of 1st, 2nd, and 3rd Order Methods", size=20)

    subfigs = fig.subfigures(1, 3)

    axs1 = subfigs[0].subplots(3, 1)
    subfigs[0].suptitle("FOG", size=18)

    x1 = data[0][:,0]
    dens1 = data[0][:,1]
    vel1 = data[0][:,2]
    pres1 = data[0][:,3]

    axs1[0].plot(x1, dens1, "--o", markerfacecolor="none", markersize=5, color="tab:red")
    axs1[0].set_xlabel("X", size = 16)
    axs1[0].set_ylabel(r"$\rho$", size=16)

    axs1[1].plot(x1, vel1, "--o", markerfacecolor="none", markersize=5, color="tab:blue")
    axs1[1].set_xlabel("X", size = 16)
    axs1[1].set_ylabel("U", size=16)

    axs1[2].plot(x1, pres1, "--o", markerfacecolor="none", markersize=5, color="tab:green")
    axs1[2].set_xlabel("X", size = 16)
    axs1[2].set_ylabel("P", size=16)


    axs2 = subfigs[1].subplots(3, 1)
    subfigs[1].suptitle("PLM + Minmod", size = 18)

    x2 = data[1][:,0]
    dens2 = data[1][:,1]
    vel2 = data[1][:,2]
    pres2 = data[1][:,3]

    axs2[0].plot(x2, dens2, "--o", markerfacecolor="none", markersize=5, color="tab:red")
    axs2[0].set_xlabel("X", size = 16)
    axs2[0].set_ylabel(r"$\rho$", size=16)

    axs2[1].plot(x2, vel2, "--o", markerfacecolor="none", markersize=5, color="tab:blue")
    axs2[1].set_xlabel("X", size = 16)
    axs2[1].set_ylabel("U", size=16)

    axs2[2].plot(x2, pres2, "--o", markerfacecolor="none", markersize=5, color="tab:green")
    axs2[2].set_xlabel("X", size = 16)
    axs2[2].set_ylabel("P", size=16)

    axs3 = subfigs[2].subplots(3, 1)
    subfigs[2].suptitle("PPM + Minmod", size=18)

    x3 = data[2][:,0]
    dens3 = data[2][:,1]
    vel3 = data[2][:,2]
    pres3 = data[2][:,3]
    
    axs3[0].plot(x3, dens3, "--o", markerfacecolor="none", markersize=5, color="tab:red")
    axs3[0].set_xlabel("X", size = 16)
    axs3[0].set_ylabel(r"$\rho$", size=16)

    axs3[1].plot(x3, vel3, "--o", markerfacecolor="none", markersize=5, color="tab:blue")
    axs3[1].set_xlabel("X", size = 16)
    axs3[1].set_ylabel("U", size=16)

    axs3[2].plot(x3, pres3, "--o", markerfacecolor="none", markersize=5, color="tab:green")
    axs3[2].set_xlabel("X", size = 16)
    axs3[2].set_ylabel("P", size=16)

    plt.savefig("./betterPlots/sodAllMethods.png", dpi=250)
    plt.close()
    

    # Second Plot
    fig, axs = plt.subplots(1, 3, constrained_layout=True, figsize=(16,8))
    fig.suptitle("Zoomed in on Pressure Discontinuity", size=20)

    titles = [
        "FOG",
        "PLM + Minmod",
        "PPM + Minmod",
    ]

    colors = [
        "tab:red",
        "tab:blue",
        "tab:green"
    ]

    for i, title in enumerate(titles):
        currData = data[i]
        x = currData[:, 0]
        dens = currData[:, 1]
        axs[i].plot(x,dens, "--o", markerfacecolor="none", markersize=5, color="tab:red")
        axs[i].set_xlim(0.4,0.9)
        axs[i].set_ylim(0.0,0.7)
        axs[i].set_xlabel("X", size=16)
        axs[i].set_ylabel(r"$\rho$", size=16)
        axs[i].set_title(title, size=20)


    plt.savefig("./betterPlots/sodDensZoom.png")
    plt.close()


    # Third Plot
    fig = plt.figure(constrained_layout=True, figsize=(16,8))
    fig.suptitle(
        "Comparison of HLL and Roe Solvers", size=20)

    subfigs = fig.subfigures(1, 2)

    primVarSyms = [
        r"$\rho$",
        "U",
        "P"
    ]
    titles = [
        "HLL",
        "Roe"
    ]
    
    colors = [
        "tab:red",
        "tab:blue",
        "tab:green"
    ]
    for j, (subfig,title) in enumerate(zip(subfigs, titles)):
        axs = subfig.subplots(3,1)
        subfig.suptitle(title, size=18)
        for i, (ax,sym, c) in enumerate(zip(axs.flatten(), primVarSyms, colors)):
            x = data[3 + j][:, 0]
            currVar = data[3 + j][:, i + 1]
            ax.plot(x, currVar, "--o", markerfacecolor="none", markersize=5, color=c)
            ax.set_xlabel("X", size=16)
            ax.set_ylabel(sym, size=16)



    plt.savefig("./betterPlots/rarefaction.png", dpi=250)    
    plt.close()


    # 4th Plot
    fig, axs = plt.subplots(1,4, constrained_layout=True, figsize=(16,8))
    fig.suptitle("Comparison of Different Slope Limiters", size=20)

    titles = [
        "Minmod",
        "MC",
        "VanLeer",
        "FOG"
    ]
    
    colors = [
        "tab:red",
        "tab:blue",
        "tab:green",
        "tab:purple"
    ]

    for i, (ax,c, title) in enumerate(zip(axs.flatten(), colors, titles)):
        x = data[5 + i][:, 0]
        dens = data[5 + i][:, 1]
        ax.plot(x, dens, "--o", markerfacecolor="none", markersize=5, color=c)
        ax.set_title(title, size=18)
        ax.set_xlabel("X", size=16)
        ax.set_ylabel(r"$\rho$", size=16)
    
    plt.savefig("./betterPlots/slopeLimiters.png", dpi=250)
    plt.close()

    # 5th Plot
    fig, axs = plt.subplots(3, 1, constrained_layout=True, figsize=(16,16) )
    fig.suptitle("Comparison of Grid Resolutions (FOG)", size=20)

    titles = [
        r"$N = 32$",
        r"$N = 64$",
        r"$N = 128$",
    ]
    
    colors = [
        "tab:red",
        "tab:blue",
        "tab:green",
    ]

    for i, (ax, c, title) in enumerate(zip(axs, colors, titles)):
        x = data[9 + i][:, 0]
        dens = data[9 + i][:, 1]
        ax.plot(x, dens, "--o", markerfacecolor="none", markersize=5, color=c)
        ax.set_title(title, size=18)
        ax.set_xlabel("X", size=16)
        ax.set_ylabel(r"$\rho$", size=16)

    plt.savefig("./betterPlots/PLMGridSizes.png", dpi=250)
    plt.close()

    # 6th Plot
    fig, axs = plt.subplots(3, 1, constrained_layout=True, figsize=(16,16) )
    fig.suptitle("Comparison of Grid Resolutions (FOG)", size=20)

    titles = [
        r"$N = 32$",
        r"$N = 64$",
        r"$N = 128$",
    ]
    
    colors = [
        "tab:red",
        "tab:blue",
        "tab:green",
    ]

    for i, (ax, c, title) in enumerate(zip(axs, colors, titles)):
        x = data[12 + i][:, 0]
        dens = data[12 + i][:, 1]
        ax.plot(x, dens, "--o", markerfacecolor="none", markersize=5, color=c)
        ax.set_title(title, size=18)
        ax.set_xlabel("X", size=16)
        ax.set_ylabel(r"$\rho$", size=16)

    plt.savefig("./betterPlots/FOGGridSizes.png", dpi=250)
    plt.close()

    # 7th Plot
    fig, axs = plt.subplots(3, 2, constrained_layout=True, figsize=(16,16) )
    fig.suptitle("Comparison of CFL Numbers", size=20)

    titles = [
        "CFL = 0.20$",
        "CFL = 0.40$",
        "CFL = 0.60",
        "CFL = 0.80",
        "CFL = 1.00",
        "CFL = 1.40",
    ]
    
    colors = [
        "tab:red",
        "tab:blue",
        "tab:green",
        "tab:purple",
        "tab:pink",
        "tab:olive"
    ]

    for i, (ax, c, title) in enumerate(zip(axs.flatten(), colors, titles)):
        x = data[15 + i][:, 0]
        dens = data[15 + i][:, 1]
        ax.plot(x, dens, "--o", markerfacecolor="none", markersize=5, color=c)
        ax.set_title(title, size=18)
        ax.set_xlabel("X", size=16)
        ax.set_ylabel(r"$\rho$", size=16)

    plt.savefig("./betterPlots/CFLComparison.png", dpi=250)
    plt.close()


