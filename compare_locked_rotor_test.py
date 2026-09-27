"""Compare Locked_rotor_test.py results with the measured Prius locked-rotor
test. Test values are read off the published plots ('Peak Locked Rotor
Torque' vs DC current, and 'Prius Locked Rotor Torque' vs rotor position),
so they carry about +/-3 N*m (peaks) to +/-5 N*m (curves) of reading error.

Figures:
  1. Peak torque vs current -- test vs simulation
  2. Torque vs electrical angle, all currents -- test (solid) vs simulation
     (dashed), same color per current
  3. One panel per current, test vs simulation
"""
import csv

import matplotlib.pyplot as plt

SIM_FILE = "torque_vs_current.csv"
SIM_KT_FILE = "torque_vs_current_kt.csv"
PEAK_PLOT_FILE = "compare_locked_rotor_peak.png"
CURVES_PLOT_FILE = "compare_locked_rotor_curves.png"
PANELS_PLOT_FILE = "compare_locked_rotor_panels.png"

# ---------------------------------------
# Test data
# ---------------------------------------
# Measured peak locked-rotor torque, DC current [A] -> torque [N*m].
TEST_PEAK = {50: 74, 75: 117, 100: 158, 125: 198, 150: 229, 200: 286, 250: 338}

# Torque [N*m] vs rotor position, 0..176 electrical deg in 8 deg steps.
TEST_ANGLES = list(range(0, 177, 8))
TEST_CURVES = {
    50: [5, -3, 0, 0, 5, 8, 13, 20, 25, 33, 45, 48, 56, 63, 72, 74, 70, 58,
         50, 46, 38, 20, 3],
    75: [0, 3, -3, -5, -3, 5, 10, 20, 37, 48, 63, 73, 82, 92, 112, 117, 115,
         102, 88, 80, 65, 35, 8],
    100: [3, -3, -8, -10, -5, 3, 8, 22, 42, 62, 82, 95, 108, 125, 150, 155,
          158, 148, 135, 120, 97, 58, 18],
    125: [3, -3, -15, -15, -12, 0, 5, 22, 50, 75, 97, 115, 132, 150, 180, 190,
          198, 192, 192, 165, 135, 80, 25],
    150: [3, -5, -25, -25, -20, -5, 5, 22, 55, 85, 110, 130, 150, 175, 210,
          222, 229, 224, 222, 198, 160, 98, 25],
    200: [5, -8, -35, -35, -25, -10, 3, 25, 62, 100, 132, 157, 190, 215, 258,
          277, 286, 280, 286, 272, 230, 135, 60],
    250: [5, -10, -38, -45, -32, -12, 5, 28, 78, 115, 160, 180, 222, 255, 303,
          325, 337, 332, 338, 331, 287, 175, 87],
}


def load_sim_curves():
    """Simulated torque vs current phase, current -> (phases, torques)."""
    curves = {}
    with open(SIM_FILE, newline="") as f:
        for r in csv.DictReader(f):
            phases, torques = curves.setdefault(float(r["Current_A"]), ([], []))
            phases.append(float(r["Phase_deg"]))
            torques.append(float(r["Torque_Nm"]))
    return curves


def load_sim_peaks():
    """Simulated peak torque, current -> (peak phase, peak torque)."""
    with open(SIM_KT_FILE, newline="") as f:
        return {
            float(r["Current_A"]): (float(r["PeakPhase_deg"]), float(r["PeakTorque_Nm"]))
            for r in csv.DictReader(f)
        }


def print_peak_table(currents, sim_peaks):
    print(f"{'I [A]':>6} {'Test [Nm]':>10} {'Sim [Nm]':>10} {'Diff [%]':>9} "
          f"{'Test ang':>9} {'Sim ang':>8}")
    for i in currents:
        t = TEST_PEAK[i]
        test_ang = TEST_ANGLES[TEST_CURVES[i].index(max(TEST_CURVES[i]))]
        sim_ang, s = sim_peaks[i]
        print(f"{i:>6} {t:>10.1f} {s:>10.1f} {100 * (s - t) / t:>+9.1f} "
              f"{test_ang:>9} {sim_ang:>8.0f}")


def plot_peaks(currents, sim_peaks):
    plt.figure()
    plt.plot(currents, [TEST_PEAK[i] for i in currents], "o-", label="Test")
    plt.plot(currents, [sim_peaks[i][1] for i in currents], "s--", label="FEMM (sector model)")
    plt.xlabel("Peak / DC Current [A]")
    plt.ylabel("Peak Locked Rotor Torque [N*m]")
    plt.title("Peak Locked Rotor Torque: Test vs Simulation")
    plt.xlim(0, 300)
    plt.ylim(0, 400)
    plt.grid(True)
    plt.legend()
    plt.savefig(PEAK_PLOT_FILE)


def plot_curves(currents, sim_curves):
    plt.figure(figsize=(9, 6))
    for k, i in enumerate(currents):
        color = f"C{k}"
        plt.plot(TEST_ANGLES, TEST_CURVES[i], "o-", color=color, markersize=3,
                 label=f"{i} A test")
        phases, torques = sim_curves[i]
        plt.plot(phases, torques, "--", color=color, label=f"{i} A sim")
    plt.xlabel("Rotor Position / Current Angle [electrical deg]")
    plt.ylabel("Torque [N*m]")
    plt.title("Prius Locked Rotor Torque: Test (solid) vs Simulation (dashed)")
    plt.xlim(0, 180)
    plt.ylim(-50, 400)
    plt.grid(True)
    plt.legend(fontsize="small", ncol=2, loc="upper left")
    plt.savefig(CURVES_PLOT_FILE)


def plot_panels(currents, sim_curves):
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.5), sharex=True)
    for ax, i in zip(axes.flat, currents):
        phases, torques = sim_curves[i]
        ax.plot(TEST_ANGLES, TEST_CURVES[i], "o-", markersize=3, label="Test")
        ax.plot(phases, torques, "s--", markersize=3, label="Sim")
        ax.set_title(f"{i} A")
        ax.set_xlim(0, 180)
        ax.grid(True)
    for ax in axes.flat[len(currents):]:
        ax.axis("off")
    axes.flat[0].legend()
    for ax in axes[1]:
        ax.set_xlabel("Electrical deg")
    for ax in axes[:, 0]:
        ax.set_ylabel("Torque [N*m]")
    fig.suptitle("Locked Rotor Torque per Current: Test vs Simulation")
    fig.tight_layout()
    fig.savefig(PANELS_PLOT_FILE)


def main():
    sim_curves = load_sim_curves()
    sim_peaks = load_sim_peaks()
    currents = [i for i in TEST_PEAK if i in sim_peaks and i in sim_curves]

    print_peak_table(currents, sim_peaks)
    plot_peaks(currents, sim_peaks)
    plot_curves(currents, sim_curves)
    plot_panels(currents, sim_curves)
    plt.show()


if __name__ == "__main__":
    main()
