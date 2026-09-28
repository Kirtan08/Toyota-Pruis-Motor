"""dq flux linkage lookup tables, and apparent vs incremental inductance.

FEMM sweep over a grid of (Id, Iq) with the magnets on: at each point the
phase flux linkages are transformed to psi_d(Id, Iq), psi_q(Id, Iq) -- the
flux maps / lookup tables a flux-based machine model or controller uses.
dq convention as in dq_inductance.py: at the locked-rotor position
(MaxTorqueInitialAngle) current angle 0 is the rotor d-axis, so
Id = I*cos(beta), Iq = I*sin(beta), and the MTPA angles of ironLoss_P2.py
give Id < 0.

Rotor position: ROTOR_POSITIONS > 1 averages the flux linkages over one slot
pitch, rotating the rotor and the current vector together (same Id, Iq),
to remove the slotting ripple from the tables.

From the tables, two kinds of inductance:
  apparent (secant):      Ld = (psi_d - psi_m) / Id,   Lq = psi_q / Iq
                          psi_m = psi_d(0, 0), the magnet flux linkage
  incremental (tangent):  Ldd = d psi_d / d Id,  Lqq = d psi_q / d Iq,
                          Ldq = d psi_d / d Iq,  Lqd = d psi_q / d Id
Apparent inductances give the flux itself (psi_d = psi_m + Ld*Id), so they
are the right ones for torque and steady-state voltage. Incremental ones
give how the flux changes with current (d psi = Ldd*dId + Ldq*dIq), so they
are the right ones for di/dt, current ripple and current-controller gains.
The comparison shows the error of using one in place of the other along
the MTPA trajectory.

iq symmetry: psi_d is even and psi_q odd in Iq (the rotor is symmetric
about the d-axis), so only Iq >= 0 is solved and the mirror is used for the
derivatives at Iq = 0.

Usage:
    python flux_map.py          # run the FEMM sweep, then plot
    python flux_map.py --plot   # replot the saved CSV only
"""
import csv
import math
import sys
import time as walltime

import femm
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

import config
import simulation
from ironLoss_P2 import MTPA_ANGLES

MODEL_FILE = "ToyotaPrius_FluxMap.FEM"
RESULTS_FILE = "flux_map.csv"
PLOT_PREFIX = "flux_map_"
RUNTIME_LOG = "flux_map_runtime.log"

# Grid, A peak: 50 A steps up to 600 A. Id > 0 only for the derivative at
# Id = 0. 14 x 13 = 182 solves (~1.8 h at ~36 s per solve). For finer
# resolution at low current, where the inductances change fastest, use 25 A
# steps below 100 A (255 solves):
#   ID_VALUES = np.concatenate([np.arange(-600, -99, 50), np.arange(-75, 51, 25)])
#   IQ_VALUES = np.concatenate([np.arange(0, 100, 25), np.arange(100, 601, 50)])
ID_VALUES = np.arange(-600, 51, 50)
IQ_VALUES = np.arange(0, 601, 50)
ROTOR_POSITIONS = 1                    # rotor positions averaged over one
                                       # slot pitch (1 = locked-rotor only)

N_SECTORS = round(360 / config.SectorAngle)   # 8 -- one pole per sector
POLE_PAIRS = config.Npoles // 2
PHASE_SHIFTS = (0, 120, 240)                  # deg, phases A, B, C
SPEED_RPM = 1000                              # for the voltage comparison


def phase_currents(i_d, i_q, theta_e):
    """Phase currents for (Id, Iq) with the d-axis at electrical angle
    theta_e -- I*sin(beta + theta_e + shift) written in Id, Iq."""
    return [i_d * simulation.sind(theta_e + s) + i_q * simulation.cosd(theta_e + s)
            for s in PHASE_SHIFTS]


def to_dq(abc, theta_e):
    """abc -> dq (amplitude invariant), inverse of phase_currents."""
    d = 2 / 3 * sum(x * simulation.sind(theta_e + s) for x, s in zip(abc, PHASE_SHIFTS))
    q = 2 / 3 * sum(x * simulation.cosd(theta_e + s) for x, s in zip(abc, PHASE_SHIFTS))
    return d, q


def solve(i_d, i_q):
    """psi_d, psi_q [Wb] (full machine) and torque [Nm], averaged over the
    rotor positions."""
    psi_d = psi_q = torque = 0.0
    for k in range(ROTOR_POSITIONS):
        rotation = k * config.ToothPitch / ROTOR_POSITIONS      # mech deg
        theta_e = POLE_PAIRS * rotation
        femm.mi_modifyboundprop(simulation.SLIDING_BAND_NAME, 10,
                                config.MaxTorqueInitialAngle + rotation)
        for name, i in zip("ABC", phase_currents(i_d, i_q, theta_e)):
            femm.mi_modifycircprop(name, 1, i)
        femm.mi_saveas(MODEL_FILE)
        femm.mi_createmesh()
        femm.mi_analyze(0)
        femm.mi_loadsolution()
        psi = [N_SECTORS * femm.mo_getcircuitproperties(n)[2] for n in "ABC"]
        d, q = to_dq(psi, theta_e)
        psi_d += d / ROTOR_POSITIONS
        psi_q += q / ROTOR_POSITIONS
        torque += femm.mo_gapintegral(simulation.SLIDING_BAND_NAME, 0) / ROTOR_POSITIONS
        femm.mo_close()
    return psi_d, psi_q, torque


def run_sweep():
    femm.openfemm()
    simulation.build_model()
    simulation.pause("Geometry, materials, and boundary conditions complete.")

    rows = []
    total = len(ID_VALUES) * len(IQ_VALUES)
    for i_d in ID_VALUES:
        for i_q in IQ_VALUES:
            psi_d, psi_q, torque = solve(float(i_d), float(i_q))
            rows.append({"Id_A": float(i_d), "Iq_A": float(i_q),
                         "PsiD_Wb": psi_d, "PsiQ_Wb": psi_q, "Torque_Nm": torque})
            print(f"[{len(rows)}/{total}] Id = {i_d:5.0f} A, Iq = {i_q:5.0f} A :: "
                  f"psi_d = {psi_d:.4f} Wb, psi_q = {psi_q:.4f} Wb, T = {torque:.1f} Nm")
    femm.closefemm()
    return rows


def save_results(rows):
    with open(RESULTS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_results():
    with open(RESULTS_FILE, newline="") as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


# ---------------------------------------
# Lookup tables and inductances
# ---------------------------------------
def tables(rows):
    """Id, Iq axes and psi_d, psi_q, torque tables indexed [id, iq]."""
    id_vals = np.array(sorted({r["Id_A"] for r in rows}))
    iq_vals = np.array(sorted({r["Iq_A"] for r in rows}))
    shape = (len(id_vals), len(iq_vals))
    psi_d, psi_q, torque = np.zeros(shape), np.zeros(shape), np.zeros(shape)
    for r in rows:
        i = np.searchsorted(id_vals, r["Id_A"])
        j = np.searchsorted(iq_vals, r["Iq_A"])
        psi_d[i, j], psi_q[i, j], torque[i, j] = r["PsiD_Wb"], r["PsiQ_Wb"], r["Torque_Nm"]
    return id_vals, iq_vals, psi_d, psi_q, torque


def inductances(id_vals, iq_vals, psi_d, psi_q):
    """Apparent and incremental inductances [H] on the table grid."""
    i0 = np.searchsorted(id_vals, 0)
    j0 = np.searchsorted(iq_vals, 0)
    psi_m = psi_d[i0, j0]

    id_grid, iq_grid = np.meshgrid(id_vals, iq_vals, indexing="ij")
    with np.errstate(invalid="ignore", divide="ignore"):
        ld_app = np.where(id_grid != 0, (psi_d - psi_m) / id_grid, np.nan)
        lq_app = np.where(iq_grid != 0, psi_q / iq_grid, np.nan)

    # Mirror to Iq < 0 (psi_d even, psi_q odd) for central differences at Iq = 0.
    iq_full = np.concatenate([-iq_vals[:0:-1], iq_vals])
    psi_d_full = np.concatenate([psi_d[:, :0:-1], psi_d], axis=1)
    psi_q_full = np.concatenate([-psi_q[:, :0:-1], psi_q], axis=1)
    keep = slice(len(iq_vals) - 1, None)
    ldd = np.gradient(psi_d_full, id_vals, axis=0)[:, keep]
    ldq = np.gradient(psi_d_full, iq_full, axis=1)[:, keep]
    lqd = np.gradient(psi_q_full, id_vals, axis=0)[:, keep]
    lqq = np.gradient(psi_q_full, iq_full, axis=1)[:, keep]
    # Limits on the axes: psi_q is odd in Iq, so psi_q/Iq -> d psi_q/d Iq at
    # Iq = 0 (exact). (psi_d - psi_m)/Id -> d psi_d/d Id only at Iq = 0; for
    # Iq > 0 cross-saturation makes it diverge, so it stays NaN there.
    lq_app[:, j0] = lqq[:, j0]
    ld_app[i0, j0] = ldd[i0, j0]
    return psi_m, {"Ld_app": ld_app, "Lq_app": lq_app,
                   "Ldd": ldd, "Lqq": lqq, "Ldq": ldq, "Lqd": lqd}


def dq_torque(psi_d, psi_q, i_d, i_q):
    return 1.5 * POLE_PAIRS * (psi_d * i_q - psi_q * i_d)


def mtpa_comparison(id_vals, iq_vals, psi_d, psi_q, torque, psi_m, ind):
    """Along the MTPA points of ironLoss_P2.py: inductances, and torque and
    1000 RPM voltage from the flux map vs from the incremental inductances
    used as if they were apparent (psi_d = psi_m + Ldd*Id, psi_q = Lqq*Iq)."""
    # Apparent L is NaN on the Id = 0 / Iq = 0 lines of the grid, so it is
    # not interpolated: it's recomputed below from the interpolated flux.
    interp = {name: RegularGridInterpolator((id_vals, iq_vals), table)
              for name, table in (("psi_d", psi_d), ("psi_q", psi_q), ("T", torque),
                                  *((k, t) for k, t in ind.items() if not k.endswith("_app")))}
    omega = POLE_PAIRS * SPEED_RPM * 2 * math.pi / 60
    out = []
    for current_amp, beta in MTPA_ANGLES.items():
        if current_amp == 0:
            continue
        i_d = current_amp * simulation.cosd(beta)
        i_q = current_amp * simulation.sind(beta)
        if not (id_vals[0] <= i_d <= id_vals[-1] and iq_vals[0] <= i_q <= iq_vals[-1]):
            continue
        v = {name: float(f((i_d, i_q))) for name, f in interp.items()}
        v["Ld_app"] = (v["psi_d"] - psi_m) / i_d
        v["Lq_app"] = v["psi_q"] / i_q
        pd_inc = psi_m + v["Ldd"] * i_d
        pq_inc = v["Lqq"] * i_q
        out.append({
            "I": current_amp, "beta": beta, "Id": i_d, "Iq": i_q, **v,
            "T_map": dq_torque(v["psi_d"], v["psi_q"], i_d, i_q),
            "T_inc": dq_torque(pd_inc, pq_inc, i_d, i_q),
            "V_map": omega * math.hypot(v["psi_d"], v["psi_q"]),
            "V_inc": omega * math.hypot(pd_inc, pq_inc),
        })
    return out


# ---------------------------------------
# Output
# ---------------------------------------
def print_summary(id_vals, iq_vals, psi_d, psi_q, torque, psi_m, ind, mtpa):
    i0 = np.searchsorted(id_vals, 0)
    print(f"Magnet flux linkage psi_m = psi_d(0, 0) = {psi_m:.4f} Wb, "
          f"psi_q(0, 0) = {psi_q[i0, 0]:.2e} Wb (should be ~0)")
    id_grid, iq_grid = np.meshgrid(id_vals, iq_vals, indexing="ij")
    t_dq = dq_torque(psi_d, psi_q, id_grid, iq_grid)
    err = np.abs(t_dq - torque)
    print(f"Torque check, 1.5*p*(psi_d*Iq - psi_q*Id) vs FEMM gap integral: "
          f"max |diff| {err.max():.1f} Nm (mean {err.mean():.1f} Nm)")
    recip = np.abs(ind["Ldq"] - ind["Lqd"])
    print(f"Reciprocity Ldq = Lqd: max |diff| {recip.max() * 1e3:.3f} mH, "
          f"max |Ldq| {np.abs(ind['Ldq']).max() * 1e3:.3f} mH")

    print(f"\nAlong MTPA (interpolated), inductances [mH]:")
    print(f"{'I [A]':>6} {'beta':>5} {'Id':>7} {'Iq':>7} {'Ld app':>7} {'Ldd':>7} "
          f"{'Lq app':>7} {'Lqq':>7} {'Ldq':>7} {'Lapp/Linc d':>11} {'q':>5}")
    for m in mtpa:
        print(f"{m['I']:>6.0f} {m['beta']:>5.0f} {m['Id']:>7.1f} {m['Iq']:>7.1f} "
              f"{m['Ld_app'] * 1e3:>7.3f} {m['Ldd'] * 1e3:>7.3f} "
              f"{m['Lq_app'] * 1e3:>7.3f} {m['Lqq'] * 1e3:>7.3f} {m['Ldq'] * 1e3:>7.3f} "
              f"{m['Ld_app'] / m['Ldd']:>11.2f} {m['Lq_app'] / m['Lqq']:>5.2f}")

    print(f"\nAlong MTPA, torque and voltage at {SPEED_RPM} RPM: flux map (= apparent L) "
          f"vs incremental L used as apparent:")
    print(f"{'I [A]':>6} {'T FEMM':>7} {'T map':>7} {'T inc':>7} {'err':>6} "
          f"{'V map':>7} {'V inc':>7} {'err':>6}")
    for m in mtpa:
        print(f"{m['I']:>6.0f} {m['T']:>7.1f} {m['T_map']:>7.1f} {m['T_inc']:>7.1f} "
              f"{100 * (m['T_inc'] / m['T_map'] - 1):>+5.0f}% "
              f"{m['V_map']:>7.1f} {m['V_inc']:>7.1f} "
              f"{100 * (m['V_inc'] / m['V_map'] - 1):>+5.0f}%")
    print("Using apparent L for di/dt instead: the current changes Lapp/Linc times "
          "slower than predicted (ratio columns above).")


def lut_tables(psi_d, psi_q, torque, ind):
    """name -> (table indexed [id, iq], unit scale, unit label)."""
    return {
        "psi_d": (psi_d, 1, "Wb"), "psi_q": (psi_q, 1, "Wb"), "torque": (torque, 1, "Nm"),
        "Ld_apparent": (ind["Ld_app"], 1e3, "mH"), "Lq_apparent": (ind["Lq_app"], 1e3, "mH"),
        "Ld_incremental": (ind["Ldd"], 1e3, "mH"), "Lq_incremental": (ind["Lqq"], 1e3, "mH"),
        "Ldq_incremental": (ind["Ldq"], 1e3, "mH"), "Lqd_incremental": (ind["Lqd"], 1e3, "mH"),
    }


def save_lut_tables(id_vals, iq_vals, luts):
    """One matrix CSV per table, ready for a lookup-table block: first row
    the Iq breakpoints, first column the Id breakpoints."""
    for name, (table, scale, unit) in luts.items():
        with open(f"{PLOT_PREFIX}lut_{name}.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([f"{name} [{unit}]: Id [A] down / Iq [A] across"]
                            + [f"{iq:g}" for iq in iq_vals])
            for i_d, row in zip(id_vals, table):
                writer.writerow([f"{i_d:g}"] + ["" if np.isnan(v) else f"{v * scale:.6g}"
                                                for v in row])


def plot_surfaces(id_vals, iq_vals, luts):
    """3D surface of each table over Id <= 0 (the motoring region). Apparent
    Ld also leaves out Id = 0, where it is undefined for Iq > 0."""
    for name, (table, scale, unit) in luts.items():
        cols = id_vals < 0 if name == "Ld_apparent" else id_vals <= 0
        id_grid, iq_grid = np.meshgrid(id_vals[cols], iq_vals, indexing="ij")
        z = table[cols] * scale
        fig = plt.figure(figsize=(10, 6.5))
        ax = fig.add_subplot(projection="3d")
        surf = ax.plot_surface(id_grid, iq_grid, z, cmap="viridis", edgecolor="k",
                               linewidth=0.3, antialiased=True)
        fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.1, label=f"{name} ({unit})")
        ax.set_xlabel("Id (A)")
        ax.set_ylabel("Iq (A)")
        ax.set_zlabel(f"{name} ({unit})")
        ax.set_title(name.replace("_", " "))
        ax.view_init(elev=25, azim=-50)
        fig.tight_layout()
        fig.savefig(f"{PLOT_PREFIX}surface_{name}.png")
        plt.close(fig)   # 9 surfaces: saved, not all opened as windows


def contour_panel(ax, id_vals, iq_vals, table, title, label, mtpa=None):
    filled = ax.contourf(id_vals, iq_vals, table.T, levels=15, cmap="viridis")
    lines = ax.contour(id_vals, iq_vals, table.T, levels=filled.levels[::2], colors="k",
                       linewidths=0.5)
    ax.clabel(lines, fontsize="x-small", fmt="%.3g")
    plt.colorbar(filled, ax=ax, label=label)
    if mtpa:
        ax.plot([m["Id"] for m in mtpa], [m["Iq"] for m in mtpa], "wo-", markersize=3,
                label="MTPA")
        ax.legend(fontsize="small", loc="upper left")
    ax.set_xlabel("Id [A]")
    ax.set_ylabel("Iq [A]")
    ax.set_title(title)


def plot_results(id_vals, iq_vals, psi_d, psi_q, torque, psi_m, ind, mtpa):
    # 1. The lookup tables.
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    contour_panel(axes[0], id_vals, iq_vals, psi_d, "psi_d(Id, Iq)", "Wb")
    contour_panel(axes[1], id_vals, iq_vals, psi_q, "psi_q(Id, Iq)", "Wb")
    contour_panel(axes[2], id_vals, iq_vals, torque, "Torque (FEMM)", "Nm", mtpa)
    fig.tight_layout()
    fig.savefig(PLOT_PREFIX + "tables.png")

    # 2. Flux linkage curves: saturation and cross-saturation.
    fig, (ax_d, ax_q) = plt.subplots(1, 2, figsize=(12, 4.8))
    for j in range(0, len(iq_vals), 2):
        ax_d.plot(id_vals, psi_d[:, j], ".-", label=f"Iq = {iq_vals[j]:g} A")
    for i in range(0, len(id_vals), 2):
        ax_q.plot(iq_vals, psi_q[i, :], ".-", label=f"Id = {id_vals[i]:g} A")
    ax_d.set_xlabel("Id [A]")
    ax_d.set_ylabel("psi_d [Wb]")
    ax_d.set_title("psi_d vs Id")
    ax_q.set_xlabel("Iq [A]")
    ax_q.set_ylabel("psi_q [Wb]")
    ax_q.set_title("psi_q vs Iq")
    for ax in (ax_d, ax_q):
        ax.grid(True)
        ax.legend(fontsize="x-small")
    fig.tight_layout()
    fig.savefig(PLOT_PREFIX + "psi_curves.png")

    # 3. Apparent vs incremental, on slices of the grid.
    fig, (ax_d, ax_q) = plt.subplots(1, 2, figsize=(12, 4.8))
    for k, j in enumerate((0, len(iq_vals) // 2, len(iq_vals) - 1)):
        ax_d.plot(id_vals, ind["Ld_app"][:, j] * 1e3, "o-", color=f"C{k}",
                  label=f"apparent, Iq = {iq_vals[j]:g} A")
        ax_d.plot(id_vals, ind["Ldd"][:, j] * 1e3, "s--", color=f"C{k}",
                  label=f"incremental, Iq = {iq_vals[j]:g} A")
    i0 = np.searchsorted(id_vals, 0)
    for k, i in enumerate((i0, i0 // 2, 0)):
        ax_q.plot(iq_vals, ind["Lq_app"][i, :] * 1e3, "o-", color=f"C{k}",
                  label=f"apparent, Id = {id_vals[i]:g} A")
        ax_q.plot(iq_vals, ind["Lqq"][i, :] * 1e3, "s--", color=f"C{k}",
                  label=f"incremental, Id = {id_vals[i]:g} A")
    ax_d.set_xlabel("Id [A]")
    ax_d.set_ylabel("Ld [mH]")
    ax_d.set_title("d-axis: apparent (solid) vs incremental (dashed)")
    ax_q.set_xlabel("Iq [A]")
    ax_q.set_ylabel("Lq [mH]")
    ax_q.set_title("q-axis: apparent (solid) vs incremental (dashed)")
    for ax in (ax_d, ax_q):
        ax.grid(True)
        ax.legend(fontsize="x-small")
    fig.tight_layout()
    fig.savefig(PLOT_PREFIX + "apparent_vs_incremental.png")

    # 4. Cross-coupling inductances.
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    contour_panel(axes[0], id_vals, iq_vals, ind["Ldq"] * 1e3, "Ldq = d psi_d / d Iq", "mH")
    contour_panel(axes[1], id_vals, iq_vals, ind["Lqd"] * 1e3, "Lqd = d psi_q / d Id", "mH")
    fig.tight_layout()
    fig.savefig(PLOT_PREFIX + "cross_coupling.png")

    # 5. Along MTPA: inductances, and torque/voltage from each model.
    if mtpa:
        currents = [m["I"] for m in mtpa]
        fig, (ax_l, ax_t, ax_v) = plt.subplots(1, 3, figsize=(16, 4.8))
        for key, style, label in (("Ld_app", "o-", "Ld apparent"), ("Ldd", "o--", "Ldd incremental"),
                                  ("Lq_app", "s-", "Lq apparent"), ("Lqq", "s--", "Lqq incremental"),
                                  ("Ldq", "^:", "Ldq cross")):
            color = "C0" if key[1] == "d" and key != "Ldq" else "C1" if key != "Ldq" else "C2"
            ax_l.plot(currents, [m[key] * 1e3 for m in mtpa], style, color=color, label=label)
        ax_l.set_ylabel("Inductance [mH]")
        ax_l.set_title("Inductances along MTPA")
        ax_t.plot(currents, [m["T"] for m in mtpa], "k.-", label="FEMM")
        ax_t.plot(currents, [m["T_map"] for m in mtpa], "o-", label="flux map / apparent L")
        ax_t.plot(currents, [m["T_inc"] for m in mtpa], "s--", label="incremental L as apparent")
        ax_t.set_ylabel("Torque [Nm]")
        ax_t.set_title("Torque along MTPA")
        ax_v.plot(currents, [m["V_map"] for m in mtpa], "o-", label="flux map / apparent L")
        ax_v.plot(currents, [m["V_inc"] for m in mtpa], "s--", label="incremental L as apparent")
        ax_v.set_ylabel("Phase voltage peak [V]")
        ax_v.set_title(f"Voltage (omega*|psi|) at {SPEED_RPM} RPM")
        for ax in (ax_l, ax_t, ax_v):
            ax.set_xlabel("Peak current [A]")
            ax.grid(True)
            ax.legend(fontsize="small")
        fig.tight_layout()
        fig.savefig(PLOT_PREFIX + "mtpa_comparison.png")


if __name__ == "__main__":
    if "--plot" in sys.argv:
        rows = load_results()
    else:
        start_time = walltime.perf_counter()
        rows = run_sweep()
        save_results(rows)
        elapsed = walltime.perf_counter() - start_time
        log_line = (f"{walltime.strftime('%Y-%m-%d %H:%M:%S')}  "
                    f"{len(ID_VALUES)}x{len(IQ_VALUES)} grid x {ROTOR_POSITIONS} positions  "
                    f"total simulation time {elapsed:.1f} s ({elapsed / 60:.1f} min)")
        print(log_line)
        with open(RUNTIME_LOG, "a") as f:
            f.write(log_line + "\n")

    id_vals, iq_vals, psi_d, psi_q, torque = tables(rows)
    psi_m, ind = inductances(id_vals, iq_vals, psi_d, psi_q)
    mtpa = mtpa_comparison(id_vals, iq_vals, psi_d, psi_q, torque, psi_m, ind)
    print_summary(id_vals, iq_vals, psi_d, psi_q, torque, psi_m, ind, mtpa)
    luts = lut_tables(psi_d, psi_q, torque, ind)
    save_lut_tables(id_vals, iq_vals, luts)
    plot_surfaces(id_vals, iq_vals, luts)
    plot_results(id_vals, iq_vals, psi_d, psi_q, torque, psi_m, ind, mtpa)
    plt.show()
