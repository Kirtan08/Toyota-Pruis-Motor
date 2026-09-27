"""Rotating-rotor sweep: magnet operating point and torque ripple.

Based on Locked_rotor_test.py, but here the rotor turns at SPEED_RPM,
synchronized with the 3-phase currents, over one electrical cycle. For each
current in CURRENTS the current angle is set to that current's MTPA angle
(the peak-torque angle from Locked_rotor_test.py), so the machine runs at
maximum torque per amp throughout.

At every step the flux density B and field strength H are read at the
middle of one magnet and projected onto the magnetization direction, giving
the magnet's operating point on its demagnetization line. The torque from
the same steps shows the torque ripple over the cycle.

Load lines: the permeance coefficient Pc is a property of the magnet and
its magnetic circuit, taken from the no-load (0 A) operating point. Stator
current doesn't change Pc; it shifts the load line B = -Pc*mu0*(H - Ha)
sideways by its demagnetizing field Ha. The CSV's 'Pc' column is the
apparent ratio -B/(mu0*H), which equals the true Pc only at no load.

Usage:
    python magnet_operating_point.py          # run the FEMM sweep, then plot
    python magnet_operating_point.py --plot   # replot the saved CSV only
"""
import csv
import math
import os
import sys
import time as walltime

import femm
import matplotlib.pyplot as plt

import config
import materials
import simulation

MODEL_FILE = "ToyotaPrius_MagnetOperatingPoint.FEM"
RESULTS_FILE = "magnet_operating_point.csv"
TORQUE_PLOT_FILE = "magnet_operating_point_torque.png"
BH_TIME_PLOT_FILE = "magnet_operating_point_BH_vs_time.png"
BH_PLOT_FILE = "magnet_operating_point_BH.png"
MTPA_FILE = "torque_vs_current_kt.csv"   # Locked_rotor_test.py peak angles

CURRENTS = [0, 125, 250]   # A, peak phase current
# MTPA current angle (electrical deg) per current, used when MTPA_FILE has
# no entry for it. The 0 A angle doesn't matter (no current).
DEFAULT_MTPA_ANGLES = {0: 0, 125: 136, 250: 144}
# Reference study's MTPA angles, from a finer (2-4 deg) current-angle sweep
# than Locked_rotor_test.py's 8 deg steps. True: use these instead.
USE_REFERENCE_MTPA = True
REFERENCE_MTPA_ANGLES = {0: 80, 125: 132, 250: 142}

# B-H plot axes, matching the reference figure.
BH_PLOT_H_MIN = -1e6   # A/m
BH_PLOT_B_MAX = 1.4    # T

SPEED_RPM = 1000
FREQ = config.Npoles * SPEED_RPM / 120   # Hz, 66.7 Hz
PERIOD = 1 / FREQ                        # s, one electrical cycle (15 ms)
NUM_STEPS = 90                           # steps per cycle -- 1 mech deg/step,
                                         # 7.5 steps per slot pitch for the ripple
STEP_TIME = PERIOD / NUM_STEPS
ROTOR_SPEED = SPEED_RPM * 360 / 60       # mech deg/s


def load_mtpa_angles():
    """MTPA angle per current: REFERENCE_MTPA_ANGLES if USE_REFERENCE_MTPA,
    else Locked_rotor_test.py's peak angles where available and
    DEFAULT_MTPA_ANGLES otherwise."""
    if USE_REFERENCE_MTPA:
        return dict(REFERENCE_MTPA_ANGLES)
    angles = dict(DEFAULT_MTPA_ANGLES)
    if os.path.exists(MTPA_FILE):
        with open(MTPA_FILE, newline="") as f:
            for r in csv.DictReader(f):
                angles[float(r["Current_A"])] = float(r["PeakPhase_deg"])
    missing = [i for i in CURRENTS if i not in angles]
    if missing:
        raise ValueError(f"No MTPA angle for currents {missing} A")
    return angles


def phase_currents(current_amp, mtpa_angle, t):
    """3-phase currents at time t: same waveform as Locked_rotor_test.py,
    advancing at FREQ from the MTPA angle."""
    return tuple(
        current_amp * simulation.sind(360 * FREQ * t + mtpa_angle + shift)
        for shift in (0, 120, 240)
    )


def magnet_probe():
    """Middle of the first magnet and the unit vector of its magnetization.
    With a sliding band the rotor mesh is never moved -- only the band's
    angle changes -- so this point stays at the magnet's middle at every
    rotor position."""
    _, magnet_pt, _ = simulation.rotor_pole_leg_label_points(0, mirror=False)
    magdir = simulation.rotor_pole_leg_magdirection(0, mirror=False)
    return magnet_pt, (simulation.cosd(magdir), simulation.sind(magdir))


def magnet_operating_point(point, direction):
    """B [T] and H [A/m] at `point`, projected onto the magnetization
    `direction`, and the apparent ratio -B / (mu0 * H) (the true
    permeance coefficient at no load)."""
    values = femm.mo_getpointvalues(*point)
    bx, by, hx, hy = values[1], values[2], values[5], values[6]
    b = bx * direction[0] + by * direction[1]
    h = hx * direction[0] + hy * direction[1]
    pc = -b / (materials.MU0 * h) if h else math.inf
    return b, h, pc


def run_rotating_sweep():
    mtpa_angles = load_mtpa_angles()
    point, direction = magnet_probe()

    femm.openfemm()
    simulation.build_model()

    simulation.pause("Geometry, materials, and boundary conditions complete.")

    rows = []
    for current_amp in CURRENTS:
        mtpa_angle = mtpa_angles[current_amp]
        for step in range(NUM_STEPS + 1):
            t = step * STEP_TIME
            rotation = ROTOR_SPEED * t
            # Rotor starts at the locked-rotor position and turns with the
            # field, so the current stays at the MTPA angle relative to it.
            femm.mi_modifyboundprop(
                simulation.SLIDING_BAND_NAME, 10,
                config.MaxTorqueInitialAngle + rotation,
            )

            i_a, i_b, i_c = phase_currents(current_amp, mtpa_angle, t)
            femm.mi_modifycircprop("A", 1, i_a)
            femm.mi_modifycircprop("B", 1, i_b)
            femm.mi_modifycircprop("C", 1, i_c)

            femm.mi_saveas(MODEL_FILE)
            femm.mi_createmesh()
            femm.mi_analyze(0)
            femm.mi_loadsolution()

            torque = femm.mo_gapintegral(simulation.SLIDING_BAND_NAME, 0)
            b, h, pc = magnet_operating_point(point, direction)
            femm.mo_close()

            rows.append({
                "Current_A": current_amp, "MTPA_deg": mtpa_angle,
                "Time_ms": t * 1e3, "RotorRotation_deg": rotation,
                "IA_A": i_a, "IB_A": i_b, "IC_A": i_c,
                "Torque_Nm": torque, "B_T": b, "H_Am": h, "Pc": pc,
            })
            print(f"I={current_amp} A :: step {step} :: {NUM_STEPS}  "
                  f"T = {torque:.1f} N*m  B = {b:.3f} T  H = {h / 1e3:.1f} kA/m")

    femm.closefemm()
    return rows


def save_results(rows):
    with open(RESULTS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_results():
    """Rows from a previous run's RESULTS_FILE, as numbers."""
    with open(RESULTS_FILE, newline="") as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def by_current(rows):
    groups = {}
    for r in rows:
        groups.setdefault(r["Current_A"], []).append(r)
    return groups


def mean(values):
    return sum(values) / len(values)


def load_lines(rows):
    """No-load permeance coefficient Pc0 (from the 0 A mean operating point)
    and, per current, the load line's shift Ha [A/m]: the line
    B = -Pc0*mu0*(H - Ha) through that current's mean operating point.
    Returns (None, {}) without a 0 A run."""
    groups = by_current(rows)
    if 0 not in groups:
        return None, {}
    b0 = mean([r["B_T"] for r in groups[0]])
    h0 = mean([r["H_Am"] for r in groups[0]])
    pc0 = -b0 / (materials.MU0 * h0)
    shifts = {}
    for current_amp, group in groups.items():
        b = mean([r["B_T"] for r in group])
        h = mean([r["H_Am"] for r in group])
        shifts[current_amp] = h + b / (pc0 * materials.MU0)
    return pc0, shifts


def print_summary(rows):
    pc0, shifts = load_lines(rows)
    if pc0 is not None:
        print(f"No-load permeance coefficient Pc = {pc0:.2f}")
    print(f"{'I [A]':>6} {'MTPA':>5} {'T mean':>8} {'Ripple pk-pk':>13} "
          f"{'B min':>7} {'B max':>7} {'H min':>9} {'Ha':>9}")
    for current_amp, group in by_current(rows).items():
        torques = [r["Torque_Nm"] for r in group]
        t_mean = mean(torques)
        ripple = max(torques) - min(torques)
        # Ripple % is meaningless at no load (mean torque ~ 0).
        ripple_pct = f"({100 * ripple / t_mean:>3.0f}%)" if abs(t_mean) > 1 else "      "
        b = [r["B_T"] for r in group]
        h = [r["H_Am"] for r in group]
        ha = f"{shifts[current_amp] / 1e3:>7.1f} k" if current_amp in shifts else ""
        print(f"{current_amp:>6.0f} {group[0]['MTPA_deg']:>5.0f} {t_mean:>8.1f} "
              f"{ripple:>7.1f} {ripple_pct} {min(b):>7.3f} {max(b):>7.3f} "
              f"{min(h) / 1e3:>7.1f} k {ha}")


def plot_results(rows):
    groups = by_current(rows)

    # Torque vs time over one electrical cycle.
    plt.figure()
    for current_amp, group in groups.items():
        plt.plot([r["Time_ms"] for r in group], [r["Torque_Nm"] for r in group],
                 ".-", label=f"{current_amp:g} A, MTPA {group[0]['MTPA_deg']:.0f} deg")
    plt.xlabel("Time, ms")
    plt.ylabel("Torque, N*m")
    plt.title(f"Rotating Torque at MTPA ({SPEED_RPM} RPM, {FREQ:.1f} Hz)")
    plt.grid(True)
    plt.legend()
    plt.savefig(TORQUE_PLOT_FILE)

    # Magnet B and H vs time.
    fig, (ax_b, ax_h) = plt.subplots(2, 1, sharex=True)
    for current_amp, group in groups.items():
        times = [r["Time_ms"] for r in group]
        ax_b.plot(times, [r["B_T"] for r in group], ".-", label=f"{current_amp:g} A")
        ax_h.plot(times, [r["H_Am"] / 1e3 for r in group], ".-", label=f"{current_amp:g} A")
    ax_b.set_ylabel("B, T")
    ax_b.set_title("Magnet Middle Point: B and H along Magnetization")
    ax_b.grid(True)
    ax_b.legend()
    ax_h.set_xlabel("Time, ms")
    ax_h.set_ylabel("H, kA/m")
    ax_h.grid(True)
    fig.savefig(BH_TIME_PLOT_FILE)

    # Operating points on the magnet's demagnetization lines.
    plt.figure()
    mu0, mur = materials.MU0, materials.MagnetMur
    for temp, hc, color in (
        (materials.MagnetRefTemp, materials.MagnetHcRef, "b"),
        (materials.MagnetTemp, materials.MagnetHc, "r"),
    ):
        plt.plot([-hc, 0], [0, mu0 * mur * hc], "-", color=color, linewidth=2,
                 label=f"N36Z_20 at {temp} C (Br = {mu0 * mur * hc:.2f} T)")
    for current_amp, group in groups.items():
        plt.plot([r["H_Am"] for r in group], [r["B_T"] for r in group],
                 "o", markersize=4, label=f"Current = {current_amp:g} Apk")
    # Load lines: no-load Pc, shifted by each current's demagnetizing field.
    pc0, shifts = load_lines(rows)
    for current_amp, ha in shifts.items():
        h_top = ha - BH_PLOT_B_MAX / (pc0 * mu0)
        plt.plot([h_top, ha], [BH_PLOT_B_MAX, 0], "k-", linewidth=1)
    if pc0 is not None:
        b0 = mean([r["B_T"] for r in groups[0]])
        h0 = mean([r["H_Am"] for r in groups[0]])
        plt.annotate(f"Pc = {pc0:.2f}", (h0, b0), xytext=(-60, 25),
                     textcoords="offset points")
    plt.xlabel("H [A/m]")
    plt.ylabel("B [T]")
    plt.title("Magnet Pc (middle of magnet)")
    plt.xlim(BH_PLOT_H_MIN, 0)
    plt.ylim(0, BH_PLOT_B_MAX)
    plt.ticklabel_format(axis="x", style="plain")
    plt.grid(True, which="both")
    plt.minorticks_on()
    plt.legend(fontsize="small", loc="upper left")
    plt.savefig(BH_PLOT_FILE)

    plt.show()


if __name__ == "__main__":
    if "--plot" in sys.argv:
        rows = load_results()
    else:
        start_time = walltime.perf_counter()

        rows = run_rotating_sweep()
        save_results(rows)

        elapsed = walltime.perf_counter() - start_time
        print(f"Simulation run time: {elapsed:.1f} s")
    print_summary(rows)

    plot_results(rows)
