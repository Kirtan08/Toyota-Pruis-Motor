"""Iron-loss data and loss-model fit for M270-35A lamination steel.

From the datasheet loss curves (W/kg vs peak B at 50-2500 Hz), fits the
loss-separation model
    P = Kh*f*B^2 + Kc*f^2*B^2 + Ke*f^1.5*B^1.5    [W/kg]
(hysteresis + classical eddy current + excess loss), and compares the
classical eddy coefficient Kc with its theoretical value for 0.35 mm
sheet, pi^2 * d^2 / (6 * rho) per unit volume.

Figures:
  1. Loss vs B per frequency (datasheet)
  2. B-H curve: bulk and corrected for the stacking factor
  3. Loss vs frequency per B
  4. P / (f * B^2) vs frequency, with the Kh + Kc*f line
  5. Loss vs B per frequency: datasheet (solid) vs model (dashed)
"""
import math

import matplotlib.pyplot as plt
import numpy as np

MU0 = 4e-7 * math.pi   # H/m

# M270-35A
THICK = 0.35        # mm, sheet thickness
BACKLACK = 0.006    # mm, insulating coating per side
SF = THICK / (THICK + 2 * BACKLACK)   # stacking factor

B_T = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
                1.3, 1.4, 1.5, 1.6, 1.7, 1.8])
H_A_M = np.array([30.0, 39.6, 46.0, 52.0, 58.2, 65.2, 73.3, 83.1, 95.5, 112,
                  136, 178, 272, 596, 1700, 3880, 7160, 11600])
# Standard M270-35A B at 2500, 5000, 10000 A/m.
H_STANDARD = [2500, 5000, 10000]
B_STANDARD = [1.49, 1.60, 1.70]

# Datasheet loss [W/kg] per frequency, at B_T[0], B_T[1], ... (each
# frequency's curve stops at a lower B the higher the frequency).
F_HZ = [50, 100, 200, 400, 1000, 2500]
P_W_KG = {
    50: [0.03, 0.07, 0.13, 0.22, 0.31, 0.43, 0.54, 0.68, 0.83, 1.01, 1.20,
         1.42, 1.70, 2.12, 2.47, 2.80, 3.05, 3.25],
    100: [0.04, 0.16, 0.34, 0.55, 0.80, 1.08, 1.38, 1.73, 2.10, 2.51, 2.98,
          3.51, 4.15, 4.97, 5.92],
    200: [0.09, 0.37, 0.79, 1.31, 1.91, 2.61, 3.39, 4.26, 5.23, 6.30, 7.51,
          8.88, 10.5, 12.5, 14.9],
    400: [0.21, 0.92, 1.99, 3.33, 4.94, 6.84, 9.00, 11.4, 14.2, 17.3, 20.9,
          24.9, 29.5, 35.4, 41.8],
    1000: [0.99, 3.67, 7.63, 12.7, 18.9, 26.4, 35.4, 46.0, 58.4, 73.0, 90.1],
    2500: [4.10, 14.9, 30.7, 52.0, 79.1, 113, 156, 209, 274, 353],
}

# Loss-model coefficients.
KH = 0.03                           # W/kg/Hz/T^2, hysteresis
KC = (0.14 - KH) / max(F_HZ)        # W/kg/Hz^2/T^2, eddy -- line through
                                    # P/(f*B^2) = 0.14 at 2500 Hz
KE = 0.0001                         # W/kg/(Hz*T)^1.5, excess
RHO = 5.2e-7                        # Ohm*m, resistivity
# Theoretical classical eddy coefficient, pi^2*d^2/(6*rho) in W/m^3, per
# kg with the reference's 1e4 divisor (~ steel density, kg/m^3).
KC_THEORY = math.pi ** 2 * (THICK / 1000) ** 2 / 6 / RHO / 1e4


def b_color(m):
    """Distinct color for the m-th B value (18 values -- more than the
    default color cycle)."""
    return plt.cm.viridis(m / (len(B_T) - 1))


def b_of(f):
    """B values that have datasheet losses at frequency f."""
    return B_T[:len(P_W_KG[f])]


def loss_model(f, b):
    return KH * f * b ** 2 + KC * f ** 2 * b ** 2 + KE * f ** 1.5 * b ** 1.5


def setup_axes(title, xlabel, ylabel):
    plt.figure()
    plt.grid(True, which="both")
    plt.minorticks_on()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


def plot_loss_vs_b():
    """Figure 1: datasheet loss vs B per frequency."""
    setup_axes("Iron Losses for M270-35A laminate steel", "B [T]", "P [W/kg]")
    for k, f in enumerate(F_HZ):
        plt.plot(b_of(f), P_W_KG[f], linewidth=2, color=f"C{k}", label=f"{f} Hz")
    plt.axis([0, B_T.max(), 0, max(P_W_KG[F_HZ[-1]])])
    plt.legend()


def plot_bh():
    """Figure 2: bulk B-H curve and the stacking-factor corrected one."""
    setup_axes("B-H curve for M270-35A laminate steel", "H [A/m]", "B [T]")
    plt.plot(H_A_M, B_T, "r", linewidth=2, label="bulk B-H curve")
    plt.plot(H_A_M, SF * B_T + (1 - SF) * H_A_M * MU0, "m", linewidth=2,
             label=f"corrected by stacking factor {SF:.4f}")
    plt.plot(H_STANDARD, B_STANDARD, "ks-", linewidth=2, markersize=5,
             label="standard M270-35A")
    plt.axis([0, H_A_M.max(), 0, B_T.max()])
    plt.legend(loc="lower right")


def plot_loss_vs_f():
    """Figure 3: loss vs frequency, one curve per B (over the frequencies
    that have data at that B)."""
    setup_axes("Iron Losses for M270-35A laminate steel", "frequency [Hz]", "P [W/kg]")
    for m, b in enumerate(B_T):
        freqs = [f for f in F_HZ if m < len(P_W_KG[f])]
        plt.plot(freqs, [P_W_KG[f][m] for f in freqs], ".-", linewidth=2,
                 color=b_color(m), label=f"{b:.1f} T")
    plt.axis([0, max(F_HZ), 0, max(P_W_KG[F_HZ[-1]])])
    plt.legend(fontsize="x-small", ncol=2, loc="upper left")


def plot_loss_per_cycle():
    """Figure 4: P / (f * B^2) vs frequency. Hysteresis gives a constant
    Kh and classical eddy loss a line of slope Kc."""
    setup_axes("Iron Losses for M270-35A laminate steel", "frequency [Hz]",
               "P/f/B^2 [W/kg/Hz/T^2]")
    for m, b in enumerate(B_T):
        freqs = [f for f in F_HZ if m < len(P_W_KG[f])]
        plt.plot(freqs, [P_W_KG[f][m] / f / b ** 2 for f in freqs], "o",
                 markersize=4, color=b_color(m), label=f"{b:.1f} T")
    plt.plot([0, max(F_HZ)], [KH, KH + KC * max(F_HZ)], "k--", linewidth=2,
             label=f"Kh + Kc*f (Kh = {KH}, Kc = {KC:.2e})")
    plt.axis([0, max(F_HZ), 0, max(P_W_KG[F_HZ[-1]]) / max(F_HZ)])
    plt.legend(fontsize="x-small", ncol=2, loc="lower right")


def plot_model_vs_data():
    """Figure 5: datasheet (solid) vs loss model (dashed)."""
    setup_axes("Iron Losses for M270-35A laminate steel", "B [T]", "P [W/kg]")
    for k, f in enumerate(F_HZ):
        plt.plot(b_of(f), P_W_KG[f], linewidth=2, color=f"C{k}", label=f"{f} Hz")
        plt.plot(B_T, loss_model(f, B_T), "--", linewidth=2, color=f"C{k}")
    plt.legend(loc="upper left")


def print_summary():
    print(f"Stacking factor SF = {SF:.4f}")
    print(f"Kh = {KH} W/kg/Hz/T^2")
    print(f"Kc = {KC:.3e} W/kg/Hz^2/T^2 (fit), {KC_THEORY:.3e} (theory, 0.35 mm)")
    print(f"Ke = {KE} W/kg/(Hz*T)^1.5")
    print(f"{'f [Hz]':>7} {'mean |err|':>11} {'max |err|':>10}   (model vs datasheet)")
    for f in F_HZ:
        data = np.array(P_W_KG[f])
        err = 100 * np.abs(loss_model(f, b_of(f)) - data) / data
        print(f"{f:>7} {err.mean():>10.0f}% {err.max():>9.0f}%")


if __name__ == "__main__":
    print_summary()
    plot_loss_vs_b()
    plot_bh()
    plot_loss_vs_f()
    plot_loss_per_cycle()
    plot_model_vs_data()
    plt.show()
