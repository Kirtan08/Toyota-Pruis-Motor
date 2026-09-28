"""Loss and efficiency maps: copper + iron loss vs speed and torque.

Combines, per current (at its MTPA angle) and speed:
  - torque: mean torque over the cycle, from ironLoss_P2.py
  - copper loss: 3/2 * Rac * Ipeak^2, Rac from phase_resistance_calculation
    at WINDING_TEMP (skin effect included)
  - iron loss: from ironLoss_P3.py (IRON_LOSS_COLUMN)
and gives
    Pout = T * omega,  Pin = Pout + Pcopper + Piron,  Eff = Pout / Pin.

Limits: the waveforms come from 1000 RPM and only frequency is scaled, and
there is no voltage limit or field weakening. Above base speed the motor
can't actually reach these operating points at MTPA -- VOLTAGE_LIMIT marks
where the phase voltage at MTPA would reach the inverter's limit. Points
with Pin above PIN_LIMIT are blanked, as in the reference.

Usage:
    python efficiency_map.py
"""
import csv
import math

import matplotlib.pyplot as plt
import numpy as np

import config
import materials
from ironLoss_P2 import by_current, load_results
from ironLoss_P2 import SUMMARY_FILE as P2_SUMMARY_FILE
from ironLoss_P3 import RESULTS_FILE as P3_RESULTS_FILE
from phase_resistance_calculation import PhaseResistDC_20C, phase_resistance_ac

PLOT_PREFIX = "efficiency_map_"
RESULTS_FILE = "efficiency_map.csv"

WINDING_TEMP = 140                  # deg C, winding temperature for Rac
IRON_LOSS_COLUMN = "PironSplit_W"   # or "PironRef_W": tooth B over the whole core
PIN_LIMIT = 60e3                    # W, points above are blanked (reference)
POWER_ISOLINES = [12.5e3, 25e3, 37.5e3, 50e3]   # W, output power; 50 kW limit

# Voltage limit: 500 V DC link (the Prius boost converter's maximum), with
# space-vector PWM giving at most Vdc/sqrt(3) peak per phase.
VDC = 500                            # V
V_PHASE_MAX = VDC / math.sqrt(3)     # V peak


def load_inputs():
    """Currents [A peak], mean torque [Nm], RPMs, iron loss [W] (current x
    RPM) and peak phase flux linkage [Wb] per current."""
    with open(P2_SUMMARY_FILE, newline="") as f:
        summary = list(csv.DictReader(f))
    currents = [float(r["Current_A"]) for r in summary]
    torque = np.array([float(r["TorqueAvg_Nm"]) for r in summary])

    with open(P3_RESULTS_FILE, newline="") as f:
        rows = list(csv.DictReader(f))
    rpms = sorted({float(r["RPM"]) for r in rows})
    iron = np.zeros((len(currents), len(rpms)))
    for r in rows:
        iron[currents.index(float(r["Current_A"])), rpms.index(float(r["RPM"]))] = \
            float(r[IRON_LOSS_COLUMN])

    groups = by_current(load_results())
    psi_peak = np.array([max(abs(r["PsiA_Wb"]) for r in groups[i]) for i in currents])
    return np.array(currents), torque, np.array(rpms), iron, psi_peak


def loss_maps(currents, torque, rpms, iron):
    r_ac = phase_resistance_ac(rpms, WINDING_TEMP)                  # per speed
    copper = 1.5 * np.outer(currents ** 2, r_ac)                    # W
    p_out = np.outer(torque, rpms * 2 * math.pi / 60)               # W
    p_in = p_out + copper + iron
    with np.errstate(invalid="ignore", divide="ignore"):
        eff = np.where(p_in > 0, 100 * p_out / p_in, 0)
    eff[p_in > PIN_LIMIT] = np.nan
    # No torque at 0 A: blank it rather than let the contours interpolate
    # from 0% up to the first loaded current.
    eff[currents == 0] = np.nan
    return copper, p_out, p_in, eff


def base_speed(currents, psi_peak):
    """Speed [RPM] at which the MTPA phase voltage, approximated as
    omega*psi_peak + Rdc*I, reaches V_PHASE_MAX."""
    r_dc = PhaseResistDC_20C * (1 + materials.CopperTempCoef * (WINDING_TEMP - 20))
    omega_per_rpm = config.Npoles / 2 * 2 * math.pi / 60        # elec rad/s per RPM
    return (V_PHASE_MAX - r_dc * currents) / (omega_per_rpm * psi_peak)


def save_results(currents, torque, rpms, iron, copper, p_out, p_in, eff):
    with open(RESULTS_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Current_A", "Torque_Nm", "RPM", "Pcopper_W", "Piron_W",
                         "Pout_W", "Pin_W", "Eff_pct"])
        for i, current_amp in enumerate(currents):
            for j, rpm in enumerate(rpms):
                writer.writerow([f"{current_amp:g}", f"{torque[i]:.2f}", f"{rpm:g}",
                                 f"{copper[i, j]:.1f}", f"{iron[i, j]:.1f}",
                                 f"{p_out[i, j]:.1f}", f"{p_in[i, j]:.1f}",
                                 f"{eff[i, j]:.2f}"])


def print_summary(currents, torque, rpms, iron, copper, eff, n_base):
    print(f"Winding at {WINDING_TEMP} C, iron loss column {IRON_LOSS_COLUMN}")
    print(f"Max copper loss {copper.max() / 1e3:.2f} kW, max iron loss {iron.max() / 1e3:.2f} kW")
    i, j = np.unravel_index(np.nanargmax(eff), eff.shape)
    print(f"Peak efficiency {eff[i, j]:.2f}% at {rpms[j]:g} RPM, "
          f"{torque[i]:.0f} Nm ({currents[i]:g} A)")
    reachable = np.where(rpms[np.newaxis, :] <= n_base[:, np.newaxis], eff, np.nan)
    i, j = np.unravel_index(np.nanargmax(reachable), eff.shape)
    print(f"  below the voltage limit: {eff[i, j]:.2f}% at {rpms[j]:g} RPM, "
          f"{torque[i]:.0f} Nm ({currents[i]:g} A)")
    print(f"\nBase speed at MTPA (Vdc {VDC} V, {V_PHASE_MAX:.0f} V peak per phase):")
    print(f"{'I [A]':>6} {'T [Nm]':>7} {'RPM':>6}")
    for current_amp, t, n in zip(currents, torque, n_base):
        print(f"{current_amp:>6.0f} {t:>7.1f} {n:>6.0f}")


def plot_map(rpms, torque, data, title, filename, levels, n_base, p_out, cbar_label):
    plt.figure(figsize=(8, 6))
    filled = plt.contourf(rpms, torque, data, levels=levels, cmap="viridis", extend="both")
    plt.colorbar(filled, label=cbar_label)
    lines = plt.contour(rpms, torque, p_out, levels=POWER_ISOLINES, colors="w",
                        linewidths=1, linestyles="--")
    plt.clabel(lines, fmt=lambda p: f"{p / 1e3:g} kW", fontsize="small")
    plt.plot(n_base, torque, "r-", linewidth=2, label="voltage limit at MTPA")
    plt.xlim(min(rpms), max(rpms))
    plt.ylim(0, max(torque))
    plt.xlabel("Speed [RPM]")
    plt.ylabel("Torque [Nm]")
    plt.title(title)
    plt.legend(loc="upper right", fontsize="small")
    plt.savefig(PLOT_PREFIX + filename)


if __name__ == "__main__":
    currents, torque, rpms, iron, psi_peak = load_inputs()
    copper, p_out, p_in, eff = loss_maps(currents, torque, rpms, iron)
    n_base = base_speed(currents, psi_peak)
    print_summary(currents, torque, rpms, iron, copper, eff, n_base)
    save_results(currents, torque, rpms, iron, copper, p_out, p_in, eff)

    kw = lambda top, n=17: np.linspace(0, top, n)
    plot_map(rpms, torque, copper / 1e3, "Winding copper loss [kW]", "copper_loss.png",
             kw(8), n_base, p_out, "kW")
    plot_map(rpms, torque, iron / 1e3, "Iron core loss [kW]", "iron_loss.png",
             kw(1.5), n_base, p_out, "kW")
    plot_map(rpms, torque, (copper + iron) / 1e3, "Total power loss [kW]", "total_loss.png",
             kw(8), n_base, p_out, "kW")
    plot_map(rpms, torque, p_out / 1e3, "Output power [kW]", "output_power.png",
             kw(60, 13), n_base, p_out, "kW")
    plot_map(rpms, torque, p_in / 1e3, "Input power [kW]", "input_power.png",
             kw(60, 13), n_base, p_out, "kW")
    plot_map(rpms, torque, eff, f"Efficiency [%] (winding {WINDING_TEMP} C)",
             "efficiency.png", np.arange(80, 97, 1), n_base, p_out, "%")
    plt.show()
