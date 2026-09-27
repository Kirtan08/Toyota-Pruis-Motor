"""Permanent-magnet B-H (demagnetization) curves for N42 NdFeB.

Reconstructs the second-quadrant B-H curve from a few datasheet values,
first with a 3-parameter model (Br, Hc, BHmax), then with a 4-parameter
model (Br, Hci, knee point Bk, recoil mur) that captures the knee, and
shifts the 4-parameter curves to other temperatures with the Br and Hci
temperature coefficients. Units: Gauss, Oersted (1 T = 1e4 Gs, 1 kA/m =
4*pi Oe).
"""
import math

import matplotlib.pyplot as plt
import numpy as np


def plot3p_bh(br, hc, bhmax, n):
    """3-parameter B-H curve (Br [Gs], Hc [Oe], BHmax [GsOe]) as n+1 points
    from H = Hc to H = 0. A hyperbola through (Hc, 0) and (0, Br) whose
    curvature 'a' is chosen so its largest B*H product equals BHmax."""
    ratio = abs(br * hc / bhmax)
    a = 2 * math.sqrt(ratio) - ratio
    h = np.linspace(hc, 0, n + 1)
    b = br * (h - hc) / (a * h - hc)
    return h, b


def plot4p_bh(br, hc, bt, mur, n):
    """4-parameter B-H curve as n+1 points from H = Hc to H = 0.

    br  -- remanence [Gs]
    hc  -- coercivity [Oe] (Hci for an intrinsic curve)
    bt  -- B [Gs] at the last point where the curve is still linear
           (the knee)
    mur -- recoil slope: mur-1 for an intrinsic curve, mur for a normal one
           (in T and A/m: (mur-1)*mu0 or mur*mu0)

    Straight line B = Br + mur*H from H = 0 down to the knee Ht, then a
    hyperbola from the knee to (Hc, 0) that meets the line at the knee."""
    ht = (bt - br) / mur
    a0 = (br - mur * (ht - hc)) / (br - mur * ht * (ht / hc - 1))
    br0 = bt * (a0 * ht - hc) / (ht - hc)
    h = np.linspace(hc, 0, n + 1)
    b = np.where(h < ht, br0 * (h - hc) / (a0 * h - hc), br + mur * h)
    return h, b


# N42 curve at -40 C
T0 = -40              # degC
Br_T = 1.38           # Tesla
Hc_Am = -1070e3       # A/m

mu0 = 4e-7 * math.pi  # H/m
mur = abs(Br_T / Hc_Am / mu0)

Br_Gs = Br_T * 1e4
Hc_Oe = Hc_Am * mu0 * 1e4
Hci_Oe = -18.7e3
Bk_Gs = 13450         # knee point of the intrinsic curve
Hk_Oe = -17000
BHimax = abs(Hk_Oe * Bk_Gs)

BHmax = 0.25 * Br_Gs ** 2 / mur   # GsOe, for a straight-line normal curve
BHmax_MGsOe = BHmax / 1e6

# Temperature coefficients
alpha = -0.110   # %/degC, Br
beta = -0.5      # %/degC, Hci
T = [20, 60, 80, 100]   # degC


def setup_axes(title):
    plt.figure()
    plt.grid(True, which="both")
    plt.minorticks_on()
    plt.axis([-20, 0, 0, 15])
    plt.xlabel("H [kOe]")
    plt.ylabel("B [kGs]")
    plt.title(title)


def plot_3_parameter():
    """Figure 1: linear curve vs. the 3-parameter curve with BHmax = 48 MGsOe."""
    setup_axes("B-H curve, 3-parameter model")
    plt.plot(np.array([Hc_Oe, 0]) / 1000, np.array([0, Br_Gs]) / 1000, label="Linear")
    h_oe, b_gs = plot3p_bh(Br_Gs, Hc_Oe, 48e6, 100)
    plt.plot(h_oe / 1000, b_gs / 1000, label="3-parameter, BHmax = 48 MGsOe")
    plt.plot(h_oe / 1000, -h_oe / 1000, "k--", label="Load line Pc = 1")
    plt.legend()


def plot_4_parameter():
    """Figure 2: 4-parameter intrinsic and normal curves at T0, then at each
    temperature in T."""
    setup_axes("B-H curve, 4-parameter model with temperature")
    plt.plot(np.array([Hc_Oe, 0]) / 1000, np.array([0, Br_Gs]) / 1000, color="gray", label="Linear")
    plt.plot([-20, 0], [20, 0], "k:", label="Load line Pc = 1")

    # Solid = intrinsic, dashed = normal; one color per temperature.
    h_oe, bi_gs = plot4p_bh(Br_Gs, Hci_Oe, Bk_Gs, mur - 1, 100)
    for k, t in enumerate([T0] + T):
        ht_oe = h_oe * (1 + beta / 100 * (t - T0))
        bit_gs = bi_gs * (1 + alpha / 100 * (t - T0))
        bnt_gs = bit_gs + ht_oe   # normal B = intrinsic B + H (CGS)
        color = f"C{k}"
        plt.plot(ht_oe / 1000, bit_gs / 1000, "-", color=color, label=f"{t} C intrinsic")
        plt.plot(ht_oe / 1000, bnt_gs / 1000, "--", color=color, label=f"{t} C normal")
    plt.legend(fontsize="small", ncol=2)


if __name__ == "__main__":
    print(f"mur = {mur:.4f}")
    print(f"Br = {Br_Gs:.0f} Gs, Hc = {Hc_Oe:.0f} Oe, Hci = {Hci_Oe:.0f} Oe")
    print(f"BHimax = {BHimax:.4g} GsOe, BHmax = {BHmax_MGsOe:.2f} MGsOe")
    plot_3_parameter()
    plot_4_parameter()
    plt.show()
