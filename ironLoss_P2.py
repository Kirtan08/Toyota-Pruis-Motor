"""Rotating-rotor sweep for iron loss (part 2): stator B waveforms.

Rotor turns at 1000 RPM synchronized with the 3-phase currents over one
electrical cycle (90 steps), as in magnet_operating_point.py, for currents
from 0 to 250 A, each at its MTPA current angle. At every step it records:
  - phase currents, flux linkages and voltages
  - torque (sliding-band gap integral)
  - B in the middle of a stator tooth and in the middle of the back iron,
    the waveforms iron losses are calculated from
  - B and H at the middle of a magnet (along its magnetization)

Voltages: in a magnetostatic solve FEMM's circuit voltage is only the
resistive drop over the modeled stack length. The induced voltage (back
EMF) is calculated here from the flux linkage, e = d(psi)/dt.

Usage:
    python ironLoss_P2.py          # run the FEMM sweep, then plot
    python ironLoss_P2.py --plot   # replot the saved CSVs only
"""
import csv
import math
import sys
import time as walltime

import femm
import matplotlib.pyplot as plt
import numpy as np

import config
import materials
import simulation
from magnet_operating_point import (
    FREQ, NUM_STEPS, ROTOR_SPEED, SPEED_RPM, STEP_TIME,
    magnet_operating_point, magnet_probe, phase_currents,
)

MODEL_FILE = "ToyotaPrius_IronLossP2.FEM"
RESULTS_FILE = "ironLoss_P2.csv"
SUMMARY_FILE = "ironLoss_P2_kt.csv"
RUNTIME_LOG = "ironLoss_P2_runtime.log"   # one line appended per run
PLOT_PREFIX = "ironLoss_P2_"

# Current [A] -> MTPA current angle [electrical deg] (reference values).
MTPA_ANGLES = {0: 90, 50: 120, 75: 124, 100: 128, 125: 132, 150: 136, 200: 140, 250: 142}
# MTPA_ANGLES = {50: 120}

N_SECTORS = round(360 / config.SectorAngle)   # 8 -- one pole per sector, in series

# B probes: middle of the stator tooth centered at one slot pitch (inside
# the sector, off its boundary), and the back iron above it.
PROBE_ANGLE = config.ToothPitch   # deg, tooth centerline
POST_HEIGHT = (config.SlotDia / 2 + config.SlotHeight
               + config.ShoeRadius + config.ShoeHeight)
TOOTH_RADIUS = config.StatorID / 2 + POST_HEIGHT / 2
BACKIRON_RADIUS = (config.StatorID / 2 + POST_HEIGHT + config.StatorOD / 2) / 2

# ORNL locked-rotor peak torque (reference script's values).
TEST_CURRENT = [50, 75, 100, 125, 150, 200, 250]
TEST_TORQUE = [74, 118, 154, 199, 221, 286, 337]


def signed_b(radius, component):
    """|B| at (radius, PROBE_ANGLE), signed by its radial ('r') or
    tangential ('t') component -- radial for the tooth, tangential for
    the back iron, where the flux runs."""
    c, s = simulation.cosd(PROBE_ANGLE), simulation.sind(PROBE_ANGLE)
    values = femm.mo_getpointvalues(radius * c, radius * s)
    bx, by = values[1], values[2]
    along = bx * c + by * s if component == "r" else -bx * s + by * c
    return math.copysign(math.hypot(bx, by), along)


def run_sweep():
    point, direction = magnet_probe()

    femm.openfemm()
    simulation.build_model()

    simulation.pause("Geometry, materials, and boundary conditions complete.")

    rows = []
    for current_amp, mtpa_angle in MTPA_ANGLES.items():
        for step in range(NUM_STEPS + 1):
            t = step * STEP_TIME
            rotation = ROTOR_SPEED * t
            femm.mi_modifyboundprop(
                simulation.SLIDING_BAND_NAME, 10,
                config.MaxTorqueInitialAngle + rotation,
            )
            for name, i in zip("ABC", phase_currents(current_amp, mtpa_angle, t)):
                femm.mi_modifycircprop(name, 1, i)

            femm.mi_saveas(MODEL_FILE)
            femm.mi_createmesh()
            femm.mi_analyze(0)
            femm.mi_loadsolution()

            row = {"Current_A": current_amp, "MTPA_deg": mtpa_angle,
                   "Time_ms": t * 1e3, "RotorRotation_deg": rotation}
            for name in "ABC":
                i, v, psi = femm.mo_getcircuitproperties(name)[:3]
                row[f"I{name}_A"] = i
                row[f"VR{name}_V"] = N_SECTORS * v     # resistive drop, full phase
                row[f"Psi{name}_Wb"] = N_SECTORS * psi  # full phase
            row["Torque_Nm"] = femm.mo_gapintegral(simulation.SLIDING_BAND_NAME, 0)
            row["Btooth_T"] = signed_b(TOOTH_RADIUS, "r")
            row["Bbackiron_T"] = signed_b(BACKIRON_RADIUS, "t")
            row["Bmag_T"], row["Hmag_Am"], _ = magnet_operating_point(point, direction)
            femm.mo_close()

            rows.append(row)
            print(f"I={current_amp} A :: step {step} :: {NUM_STEPS}  "
                  f"T = {row['Torque_Nm']:.1f} N*m  Btooth = {row['Btooth_T']:.2f} T  "
                  f"Bbackiron = {row['Bbackiron_T']:.2f} T")

    femm.closefemm()
    add_back_emf(rows)
    add_terminal_voltage(rows)
    return rows


def add_back_emf(rows):
    """EMF = d(psi)/dt per phase, from the flux linkage over the cycle
    (central differences; the cycle is periodic, so the last step, which
    repeats the first, is dropped before differencing)."""
    for group in by_current(rows).values():
        for name in "ABC":
            psi = np.array([r[f"Psi{name}_Wb"] for r in group[:-1]])
            emf = (np.roll(psi, -1) - np.roll(psi, 1)) / (2 * STEP_TIME)
            emf = np.append(emf, emf[0])
            for r, e in zip(group, emf):
                r[f"EMF{name}_V"] = e


def add_terminal_voltage(rows):
    """Phase terminal voltage v = R*i + d(psi)/dt: FEMM's resistive drop
    (active length only, no end turns) plus the induced EMF."""
    for r in rows:
        for name in "ABC":
            r[f"VT{name}_V"] = r[f"VR{name}_V"] + r[f"EMF{name}_V"]


def by_current(rows):
    groups = {}
    for r in rows:
        groups.setdefault(r["Current_A"], []).append(r)
    return groups


def summarize(rows):
    """Per current: MTPA angle and max / mean / min torque. The mean skips
    the last step, which repeats the first."""
    summary = []
    for current_amp, group in by_current(rows).items():
        torques = [r["Torque_Nm"] for r in group]
        summary.append({
            "Current_A": current_amp, "MTPA_deg": group[0]["MTPA_deg"],
            "TorqueMax_Nm": max(torques),
            "TorqueAvg_Nm": sum(torques[:-1]) / len(torques[:-1]),
            "TorqueMin_Nm": min(torques),
            "BtoothPeak_T": max(abs(r["Btooth_T"]) for r in group),
            "BbackironPeak_T": max(abs(r["Bbackiron_T"]) for r in group),
            "EMFpeak_V": max(abs(r["EMFA_V"]) for r in group),
        })
    return summary


def save_results(rows, summary):
    for path, data in ((RESULTS_FILE, rows), (SUMMARY_FILE, summary)):
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)


def load_results():
    with open(RESULTS_FILE, newline="") as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def print_summary(summary):
    print(f"{'I [A]':>6} {'MTPA':>5} {'T max':>7} {'T avg':>7} {'T min':>7} "
          f"{'B tooth':>8} {'B back':>7} {'EMF pk':>7}")
    for s in summary:
        print(f"{s['Current_A']:>6.0f} {s['MTPA_deg']:>5.0f} {s['TorqueMax_Nm']:>7.1f} "
              f"{s['TorqueAvg_Nm']:>7.1f} {s['TorqueMin_Nm']:>7.1f} "
              f"{s['BtoothPeak_T']:>8.2f} {s['BbackironPeak_T']:>7.2f} {s['EMFpeak_V']:>7.1f}")


def plot_phase_quantity(groups, key, ylabel, title, filename):
    """One color per current; phase A solid, B dashed, C dotted."""
    plt.figure(figsize=(8, 5))
    for k, (current_amp, group) in enumerate(groups.items()):
        times = [r["Time_ms"] for r in group]
        for name, style in zip("ABC", ("-", "--", ":")):
            plt.plot(times, [r[key.format(name)] for r in group], style, color=f"C{k}",
                     label=f"{current_amp:g} A" if name == "A" else None)
    plt.xlabel("Time [ms]")
    plt.ylabel(ylabel)
    plt.title(f"{title} (A solid, B dashed, C dotted)")
    plt.grid(True)
    plt.legend(fontsize="small", ncol=2)
    plt.savefig(PLOT_PREFIX + filename)


def plot_per_current(groups, key, ylabel, title, filename):
    plt.figure(figsize=(8, 5))
    for k, (current_amp, group) in enumerate(groups.items()):
        plt.plot([r["Time_ms"] for r in group], [r[key] for r in group], ".-",
                 color=f"C{k}", label=f"{current_amp:g} A")
    plt.xlabel("Time [ms]")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.legend(fontsize="small", ncol=2)
    plt.savefig(PLOT_PREFIX + filename)


def plot_results(rows, summary):
    groups = by_current(rows)

    plot_phase_quantity(groups, "I{}_A", "Current [A]", "Phase Currents", "currents.png")
    plot_phase_quantity(groups, "EMF{}_V", "Induced Voltage [V]",
                        f"Induced Phase Voltages at {SPEED_RPM} RPM", "emf.png")
    # What the reference script plots as "Induced Phase Voltages": FEMM's
    # magnetostatic circuit voltage, i.e. only the resistive drop.
    plot_phase_quantity(groups, "VR{}_V", "Resistive Voltage [V]",
                        "FEMM Circuit Voltage (R*i)", "voltage_resistive.png")
    plot_phase_quantity(groups, "VT{}_V", "Terminal Voltage [V]",
                        f"Phase Terminal Voltages R*i + dpsi/dt at {SPEED_RPM} RPM",
                        "voltage_terminal.png")
    plot_phase_quantity(groups, "Psi{}_Wb", "Flux Linkage [Wb]",
                        "Phase Flux Linkages", "flux_linkage.png")
    plot_per_current(groups, "Torque_Nm", "Torque [Nm]", "Transient Torque", "torque.png")
    plot_per_current(groups, "Btooth_T", "B field tooth [T]",
                     "Stator Tooth B Field (radial)", "B_tooth.png")
    plot_per_current(groups, "Bbackiron_T", "B field back-iron [T]",
                     "Stator Back-Iron B Field (tangential)", "B_backiron.png")

    # Magnet operating points on the demagnetization lines.
    plt.figure()
    mu0, mur = materials.MU0, materials.MagnetMur
    for k, (current_amp, group) in enumerate(groups.items()):
        plt.plot([r["Hmag_Am"] for r in group], [r["Bmag_T"] for r in group], ".",
                 color=f"C{k}", label=f"{current_amp:g} A")
    for temp, hc, color in ((materials.MagnetRefTemp, materials.MagnetHcRef, "b"),
                            (materials.MagnetTemp, materials.MagnetHc, "r")):
        plt.plot([-hc, 0], [0, mu0 * mur * hc], "-", color=color, linewidth=2,
                 label=f"N36Z_20 {temp} C")
    plt.xlabel("H [A/m]")
    plt.ylabel("B [T]")
    plt.title("Magnet Pc")
    plt.axis([-1e6, 0, 0, 1.4])
    plt.ticklabel_format(axis="x", style="plain")
    plt.grid(True)
    plt.legend(fontsize="small", loc="upper left", ncol=2)
    plt.savefig(PLOT_PREFIX + "magnet_BH.png")

    # Torque: test vs simulation max / average / min.
    plt.figure()
    currents = [s["Current_A"] for s in summary]
    plt.plot(TEST_CURRENT, TEST_TORQUE, ".-", color="b", linewidth=2, label="Test")
    plt.plot(currents, [s["TorqueMax_Nm"] for s in summary], "--", color="r",
             linewidth=2, label="Simulation Max")
    plt.plot(currents, [s["TorqueAvg_Nm"] for s in summary], ".-", color="r",
             linewidth=2, label="Simulation Avg")
    plt.plot(currents, [s["TorqueMin_Nm"] for s in summary], ":", color="r",
             linewidth=2, label="Simulation Min")
    plt.xlabel("Current [A]")
    plt.ylabel("Torque [Nm]")
    plt.title("Torque Test vs Simulation")
    plt.axis([0, max(TEST_CURRENT), 0, 400])
    plt.grid(True)
    plt.legend(loc="upper left")
    plt.savefig(PLOT_PREFIX + "torque_vs_test.png")

    plt.show()


if __name__ == "__main__":
    if "--plot" in sys.argv:
        rows = load_results()
        if "VTA_V" not in rows[0]:   # CSVs saved before the column existed
            add_terminal_voltage(rows)
    else:
        start_time = walltime.perf_counter()
        rows = run_sweep()
        elapsed = walltime.perf_counter() - start_time
        log_line = (f"{walltime.strftime('%Y-%m-%d %H:%M:%S')}  "
                    f"currents {list(MTPA_ANGLES)} A  "
                    f"{len(MTPA_ANGLES) * (NUM_STEPS + 1)} solves  "
                    f"total simulation time {elapsed:.1f} s ({elapsed / 60:.1f} min)")
        print(log_line)
        with open(RUNTIME_LOG, "a") as f:
            f.write(log_line + "\n")
    summary = summarize(rows)
    if "--plot" not in sys.argv:
        save_results(rows, summary)
    print_summary(summary)
    plot_results(rows, summary)
