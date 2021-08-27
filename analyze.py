#!/bin/python3

import matplotlib.pyplot as plt
import load

def titrate(lambdaFileName, cutoff=0.80):
    lambda_proto   = 0
    lambda_deproto = 0
    
    for x in load.Col(lambdaFileName, 2):
        if (x > cutoff):
            lambda_deproto += 1
        if (x < 1 - cutoff):
            lambda_proto   += 1

    fraction = float(lambda_deproto) / (lambda_proto + lambda_deproto)

    return fraction

def compareLambdaFiles(namelist):
    # If you accidentally put a string instead of a list, fix it.
    if (type(namelist) == type("")):
        namelist = [namelist]

    # Define (size of) main figure.
    fig = plt.figure(figsize=(24, 10))

    # Define sub-plots.
    plt1 = fig.add_subplot(2, 4, 1)
    plt2 = fig.add_subplot(2, 4, 2)
    plt3 = fig.add_subplot(2, 4, 3)
    plt4 = fig.add_subplot(2, 4, 4)
    plt5 = fig.add_subplot(2, 4, 5)
    plt6 = fig.add_subplot(2, 4, 6)
    plt7 = fig.add_subplot(2, 4, 7)
    plt8 = fig.add_subplot(2, 4, 8)

    # Get the data and plot.
    for name in namelist:
        time        = load.Col(name, 1)
        lambda_x    = load.Col(name, 2)
        lambda_dvdl = load.Col(name, 3)
        lambda_temp = load.Col(name, 4)
        lambda_vel  = load.Col(name, 5)
        F_coulomb   = load.Col(name, 6)
        F_corr      = load.Col(name, 7)
        F_bias      = load.Col(name, 8)
        F_ph        = load.Col(name, 9)

        plt1.plot(time, lambda_x, linewidth=0.5, label="deprotonation = {:.2f}".format(titrate(name)))
        plt2.plot(time, lambda_temp, linewidth=0.5, label="mean = {:.1f} (K)".format(sum(lambda_temp)/len(lambda_temp)))
        plt3.hist(lambda_vel, density=True)
        plt4.scatter(lambda_x, lambda_dvdl, s=5)
        plt5.scatter(lambda_x, F_coulomb, s=5)
        plt6.scatter(lambda_x, F_corr, s=5)
        plt7.scatter(lambda_x, F_bias, s=5)
        plt8.scatter(lambda_x, F_ph, s=5)

    plt1.set_title("$\lambda$-coordinate vs time")
    plt1.set_xlabel("Time (ps)")
    plt1.set_ylabel("$\lambda$-coordinate")
    plt1.set_ylim(-0.1, 1.1)
    plt1.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
    plt1.legend()

    plt2.set_title("$\lambda$-temperature vs time")
    plt2.set_xlabel("Time (ps)")
    plt2.set_ylabel("$\lambda$-temperature (K)")
    plt2.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
    plt2.legend()

    plt3.set_title("$\lambda$-velocity distribution")
    plt3.set_xlabel("$\lambda$-velocity (km/s)")

    plt4.set_title("Force (dV/dl) on $\lambda$-particle")
    plt4.set_xlabel("$\lambda$-coordinate")
    plt4.set_ylabel("dV/dl")
    plt4.set_xlim(-0.1, 1.1)

    plt5.set_title("Coulomb-force on $\lambda$-particle")
    plt5.set_xlabel("$\lambda$-coordinate")
    plt5.set_ylabel("$F_{Coulomb}$")
    plt5.set_xlim(-0.1, 1.1)

    plt6.set_title("Reference-force on $\lambda$-particle")
    plt6.set_xlabel("$\lambda$-coordinate")
    plt6.set_ylabel("$F_{corr}$")
    plt6.set_xlim(-0.1, 1.1)

    plt7.set_title("Bias-force on $\lambda$-particle")
    plt7.set_xlabel("$\lambda$-coordinate")
    plt7.set_ylabel("$F_{bias}$")
    plt7.axis([-0.1, 1.1, -200, 200])

    plt8.set_title("pH-force on $\lambda$-particle")
    plt8.set_xlabel("$\lambda$-coordinate")
    plt8.set_ylabel("$F_{pH}$")
    plt8.set_xlim(-0.1, 1.1)

    # Stuff we do in all subplots we can do in a loop:
    for plot in [plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8]:
        plot.grid()

    fig.legend(namelist, loc="upper center")
    # plt.tight_layout()
    plt.show()

compareLambdaFiles(['lambda_1.dat', 'lambda_2.dat', 'lambda_3.dat'])
compareLambdaFiles(['lambda_4.dat', 'lambda_5.dat', 'lambda_6.dat'])
