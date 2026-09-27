"""Analytical dq-axis torque model of the Prius IPM motor vs the ORNL
locked-rotor test.

Torque of an interior-PM machine in the dq frame:
    T = 3/2 * pp * (Psi_PM * Iq + (Ld - Lq) * Id * Iq)
      = magnet torque  3/2 * pp * Psi_PM * I * cos(gamma)
      + reluctance torque  -3/4 * pp * (Ld - Lq) * I^2 * sin(2*gamma)
with Id = -I*sin(gamma), Iq = I*cos(gamma), gamma = current angle from the
q-axis (gamma = beta - 90, beta being the test's rotor position in
electrical degrees). Psi_PM and (Ld - Lq) are fitted per current, so they
include saturation.

Figures:
  1. Test torque vs gamma (dashed) with the total dq-model torque (black)
     and its MTPA points
  2. Test torque vs gamma with the magnet (magenta) and reluctance (red)
     torque components
  3. Torque vs Iq: test peaks vs the dq model
"""
import matplotlib.pyplot as plt
import numpy as np

POLES = 8
PP = POLES // 2   # pole pairs

beta = np.arange(0, 177, 8)   # electrical deg (rotor position in the test)
gamma = beta - 90             # current angle from the q-axis, deg

# ORNL locked-rotor test torque [N*m] at each beta, per current [A].
CURRENTS = [100, 150, 200, 250]
TEST_TORQUE = {
    250: [10, -10.1, -43.3, -46.2, -33.8, -16.4, -1.9, 31.7, 78.1, 114.2, 158.5,
          181.6, 222, 255, 302, 324, 337, 332, 339, 331, 287, 175, 87],
    200: [6, -9.6, -37.6, -35.4, -26.7, -11.7, 0, 27, 62, 99.2, 132, 157, 189.2,
          214.3, 258.6, 277, 286.1, 280.4, 286.6, 271.2, 228.2, 134.2, 58.3],
    150: [3.6, -7.6, -25.7, -24.7, -17.3, -3.5, 7.4, 25.8, 54.8, 83.6, 109.9,
          127.8, 151, 177.4, 208.1, 223, 229, 223, 221, 198.2, 159.6, 99.7, 25.3],
    100: [3.5, 0.7, -8.1, -9, -5.8, 6.1, 10.4, 22, 41, 61, 82, 95, 109, 125, 149,
          157, 158, 149, 136, 120, 98, 59, 20.4],
}

# dq-model parameters, one per current in CURRENTS (fitted to the test).
PSI_PM = [0.175, 0.155, 0.135, 0.125]           # Wb, magnet flux linkage
LD_MINUS_LQ = [-0.0025, -0.0018, -0.00136, -0.00105]   # H, Ld - Lq
LD = 0.00196   # H, unsaturated d-axis inductance (reference only)
LQ = 0.00347   # H, unsaturated q-axis inductance (reference only)

# ORNL peak torque vs current (Kt test) and its rotor position.
CURRENT_MAX = np.array([50, 75, 100, 125, 150, 200, 250])
TORQUE_MAX = np.array([74, 118, 154, 199, 221, 286, 337])
BETA_MAX = np.full(len(CURRENT_MAX), 128)
GAMMA_MAX = BETA_MAX - 90


def sind(x):
    return np.sin(np.radians(x))


def cosd(x):
    return np.cos(np.radians(x))


def dq_torque(current, psi_pm, ld_minus_lq):
    """Total, magnet and reluctance torque vs gamma for one current."""
    i_d = -current * sind(gamma)
    i_q = current * cosd(gamma)
    total = 1.5 * PP * (psi_pm * i_q + ld_minus_lq * i_d * i_q)
    magnet = 1.5 * PP * psi_pm * current * cosd(gamma)
    reluctance = 1.5 * PP * (-0.5 * ld_minus_lq * current ** 2 * sind(2 * gamma))
    return total, magnet, reluctance


def compute():
    """dq-model torques per current and the MTPA point (max total torque
    over the gamma grid) of each."""
    results = {}
    for current, psi_pm, dl in zip(CURRENTS, PSI_PM, LD_MINUS_LQ):
        total, magnet, reluctance = dq_torque(current, psi_pm, dl)
        k = int(np.argmax(total))
        results[current] = {
            "total": total, "magnet": magnet, "reluctance": reluctance,
            "mtpa_torque": total[k], "mtpa_gamma": gamma[k],
            "mtpa_id": current * sind(gamma[k]),
            "mtpa_iq": current * cosd(gamma[k]),
        }
    return results


def setup_torque_axes(title):
    plt.figure(figsize=(8, 6))
    plt.grid(True, which="both")
    plt.minorticks_on()
    plt.axis([-90, 90, -50, 350])
    plt.xlabel("Gamma - current angle [deg]")
    plt.ylabel("Torque [Nm]")
    plt.title(title)


def plot_test(colors):
    for current in CURRENTS:
        plt.plot(gamma, TEST_TORQUE[current], "--", linewidth=2,
                 color=colors[current], label=f"Test {current} A")


def plot_results(results):
    colors = {current: f"C{k}" for k, current in enumerate(CURRENTS)}

    # Figure 1: total dq torque and MTPA trajectory.
    setup_torque_axes("Torque vs gamma - together magnet torque and reluctance torque")
    plot_test(colors)
    for k, current in enumerate(CURRENTS):
        plt.plot(gamma, results[current]["total"], linewidth=2, color="k",
                 label="dq model" if k == 0 else None)
    plt.plot([results[c]["mtpa_gamma"] for c in CURRENTS],
             [results[c]["mtpa_torque"] for c in CURRENTS],
             "o-", linewidth=2, color="k", markerfacecolor="w", label="MTPA")
    plt.legend(loc="upper left")

    # Figure 2: magnet and reluctance torque components.
    setup_torque_axes("Torque vs gamma - separate magnet torque and reluctance torque")
    plot_test(colors)
    for k, current in enumerate(CURRENTS):
        plt.plot(gamma, results[current]["magnet"], linewidth=2, color="m",
                 label="Magnet torque" if k == 0 else None)
        plt.plot(gamma, results[current]["reluctance"], linewidth=2, color="r",
                 label="Reluctance torque" if k == 0 else None)
    plt.legend(loc="upper left")

    # Figure 3: torque vs Iq. The dq-model line uses the torque at gamma =
    # -2 deg (beta = 88), where Id ~ 0 and the torque is almost pure magnet
    # torque, 3/2 * pp * Psi_PM * Iq.
    k_near_q = 11
    plt.figure()
    plt.grid(True, which="both")
    plt.minorticks_on()
    plt.axis([0, 250, 0, 350])
    plt.xlabel("Current Iq [A]")
    plt.ylabel("Torque [Nm]")
    plt.title("Torque vs current - KT")
    plt.plot(np.r_[0, CURRENT_MAX * cosd(GAMMA_MAX)], np.r_[0, TORQUE_MAX],
             linewidth=2, color="r", label="KT = Tmax/Iq (test)")
    plt.plot([0] + [results[c]["mtpa_iq"] for c in CURRENTS],
             [0] + [results[c]["total"][k_near_q] for c in CURRENTS],
             linewidth=2, color="b", label="KT = 3/2*p*Psi_PM (dq model)")
    plt.legend(loc="upper left")

    plt.show()


def print_results(results):
    print(f"{'I [A]':>6} {'Psi_PM':>7} {'Ld-Lq [mH]':>11} {'MTPA T':>8} "
          f"{'gamma':>6} {'beta':>5} {'Id':>7} {'Iq':>7}")
    for current, psi_pm, dl in zip(CURRENTS, PSI_PM, LD_MINUS_LQ):
        r = results[current]
        print(f"{current:>6} {psi_pm:>7.3f} {dl * 1e3:>11.2f} {r['mtpa_torque']:>8.1f} "
              f"{r['mtpa_gamma']:>6.0f} {r['mtpa_gamma'] + 90:>5.0f} "
              f"{r['mtpa_id']:>7.1f} {r['mtpa_iq']:>7.1f}")


if __name__ == "__main__":
    results = compute()
    print_results(results)
    plot_results(results)
