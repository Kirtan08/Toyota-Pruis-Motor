"""d- and q-axis inductances vs current.

Same setup as the reference inductance study: magnets replaced by air, so
only the stator current makes flux; rotor fixed at the locked-rotor
position (MaxTorqueInitialAngle), where current angle 0 deg lines the
current up with the rotor d-axis and 90 deg with the q-axis (as in
Locked_rotor_test.py: Id = I*cos(beta), Iq = I*sin(beta)).

For each current the model is solved once on each axis. The phase flux
linkages from FEMM's circuit results are transformed to psi_d, psi_q with
the same transform the currents use, and
    Ld = psi_d / Id  (current on the d-axis)
    Lq = psi_q / Iq  (current on the q-axis)
The sector holds one pole of each phase winding; the poles are in series,
so the full-machine flux linkage is N_SECTORS times the sector's.

Also saves the airgap radial flux density along the sector for each run.

Magnetizing vs leakage inductance: the circuit flux linkage counts every
flux path, including slot and tooth-tip leakage. The magnetizing part is
the flux linkage of the airgap flux alone -- the radial B across the
airgap integrated against each phase's winding (turns) function:
    psi_m = Length * r * integral( Br(theta) * W(theta) dtheta )
Lmd, Lmq come from psi_m the same way as Ld, Lq, and the leakage inductance
is L - Lm. Valid with MAGNETS_AS_AIR = True (no magnet flux in Br).

Usage:
    python dq_inductance.py          # run the FEMM sweep, then plot
    python dq_inductance.py --plot   # replot the saved CSVs only
"""
import csv
import sys
import time as walltime

import femm
import matplotlib.pyplot as plt
import numpy as np

import config
import materials
import simulation

MODEL_FILE = "ToyotaPrius_DqInductance.FEM"
RESULTS_FILE = "dq_inductance.csv"
AIRGAP_FILE = "dq_inductance_airgap_B.csv"
INDUCTANCE_PLOT_FILE = "dq_inductance.png"
AIRGAP_PLOT_FILE = "dq_inductance_airgap_B.png"
MAGNETIZING_PLOT_FILE = "dq_inductance_magnetizing.png"

CURRENTS = [50, 75, 100, 125, 150, 200, 250]   # A, peak phase current
AXES = {"d": 0, "q": 90}                        # current angle, electrical deg

# True: magnets set to air (the reference's method), giving the inductances
# of the stator field alone. False: magnets kept, and the no-load flux
# linkage is subtracted on the d-axis (apparent inductances).
MAGNETS_AS_AIR = True

N_SECTORS = round(360 / config.SectorAngle)   # 8 -- one pole per sector
PHASE_SHIFTS = (0, 120, 240)                  # deg, phases A, B, C

# Reference values from torque_calculation_PM+RELUCTANCE.PY.
REF_LD = 1.96e-3   # H, unsaturated
REF_LQ = 3.47e-3   # H, unsaturated
REF_LD_MINUS_LQ = {100: -2.5e-3, 150: -1.8e-3, 200: -1.36e-3, 250: -1.05e-3}  # H, fitted
# Reference study's 'D-Q axes inductances' chart, read off the plot
# (about +/-0.05 mH), current [A] -> H.
REF_CHART_LD = {50: 1.85e-3, 75: 1.65e-3, 100: 1.53e-3, 125: 1.43e-3,
                150: 1.35e-3, 200: 1.22e-3, 250: 1.10e-3}
REF_CHART_LQ = {50: 4.55e-3, 75: 3.55e-3, 100: 2.87e-3, 125: 2.40e-3,
                150: 2.07e-3, 200: 1.62e-3, 250: 1.33e-3}

INCH = 0.0254   # m per inch -- the model is in inches

# Airgap B sampling: radius between the sliding band and the stator bore.
AIRGAP_RADIUS = (simulation.SLIDING_BAND_RADIUS + config.StatorID / 2) / 2
AIRGAP_POINTS = 181


def phase_currents(current_amp, beta):
    return [current_amp * simulation.sind(beta + s) for s in PHASE_SHIFTS]


def to_dq(a, b, c):
    """abc -> dq with the transform that maps the currents above to
    Id = I*cos(beta), Iq = I*sin(beta)."""
    abc = (a, b, c)
    d = 2 / 3 * sum(x * simulation.sind(s) for x, s in zip(abc, PHASE_SHIFTS))
    q = 2 / 3 * sum(x * simulation.cosd(s) for x, s in zip(abc, PHASE_SHIFTS))
    return d, q


def winding_function(phase, theta):
    """Turns function W(theta) of one phase over one pole pair (two sectors,
    the second the anti-periodic mirror of the first), zero mean. Each slot
    adds its signed turns as theta passes its center. The sign is set so a
    phase's airgap flux linkage has the same sign as FEMM's circuit flux
    linkage for that phase."""
    circuits = config.SlotCircuits * 2
    dirs = config.SlotCoilDirs + [-d for d in config.SlotCoilDirs]
    w = np.zeros_like(theta)
    for k, (circuit, direction) in enumerate(zip(circuits, dirs)):
        if circuit == phase:
            center = config.ToothPitch * (k + 0.5)
            w -= direction * config.Turns * (theta >= center)
    return w - w.mean()


def magnetizing_flux_linkages(angles, br):
    """Full-machine flux linkage [Wb-turn] of phases A, B, C from the
    airgap flux alone. `br` covers one sector; the next sector is its
    anti-periodic image (Br -> -Br)."""
    angles = np.asarray(angles, dtype=float)
    br = np.asarray(br, dtype=float)
    theta = np.concatenate([angles, angles[1:] + config.SectorAngle])
    b = np.concatenate([br, -br[1:]])
    pole_pairs = N_SECTORS // 2
    r = AIRGAP_RADIUS * INCH
    length = config.Length * INCH
    return [
        pole_pairs * length * r * np.trapz(b * winding_function(phase, theta), np.radians(theta))
        for phase in "ABC"
    ]


def magnetizing_inductances(airgap):
    """current -> (Lmd, Lmq) in H, from the saved airgap Br."""
    table = {}
    for current_amp, axis, angles, br in airgap:
        psi_d, psi_q = to_dq(*magnetizing_flux_linkages(angles, br))
        i_d, i_q = to_dq(*phase_currents(current_amp, AXES[axis]))
        lm = psi_d / i_d if axis == "d" else psi_q / i_q
        table.setdefault(current_amp, {})[axis] = lm
    return {i: (v["d"], v["q"]) for i, v in table.items()}


def phase_flux_linkages():
    """Full-machine flux linkage [Wb-turn] of phases A, B, C."""
    return [N_SECTORS * femm.mo_getcircuitproperties(name)[2] for name in "ABC"]


def airgap_br():
    """Radial flux density [T] along the airgap arc across the sector."""
    angles = np.linspace(0, config.SectorAngle, AIRGAP_POINTS)
    br = []
    for a in angles:
        x, y = AIRGAP_RADIUS * simulation.cosd(a), AIRGAP_RADIUS * simulation.sind(a)
        values = femm.mo_getpointvalues(x, y)
        br.append(values[1] * simulation.cosd(a) + values[2] * simulation.sind(a))
    return angles, br


def solve(current_amp, beta):
    """Solve with the given current on the given angle; return the phase
    flux linkages, torque and airgap B."""
    for name, i in zip("ABC", phase_currents(current_amp, beta)):
        femm.mi_modifycircprop(name, 1, i)
    femm.mi_saveas(MODEL_FILE)
    femm.mi_createmesh()
    femm.mi_analyze(0)
    femm.mi_loadsolution()
    psi = phase_flux_linkages()
    torque = femm.mo_gapintegral(simulation.SLIDING_BAND_NAME, 0)
    angles, br = airgap_br()
    femm.mo_close()
    return psi, torque, angles, br


def run_inductance_sweep():
    femm.openfemm()
    simulation.build_model()
    if MAGNETS_AS_AIR:
        # Hc = 0 and mu = 1: the magnet regions behave as air.
        femm.mi_modifymaterial(materials.MAGNET, 1, 1)
        femm.mi_modifymaterial(materials.MAGNET, 2, 1)
        femm.mi_modifymaterial(materials.MAGNET, 3, 0)

    simulation.pause("Geometry, materials, and boundary conditions complete.")

    femm.mi_modifyboundprop(
        simulation.SLIDING_BAND_NAME, 10, config.MaxTorqueInitialAngle
    )

    psi_d0 = 0.0
    if not MAGNETS_AS_AIR:
        psi0, _, _, _ = solve(0, 0)
        psi_d0, _ = to_dq(*psi0)

    rows = []
    airgap = []
    for current_amp in CURRENTS:
        for axis, beta in AXES.items():
            psi, torque, angles, br = solve(current_amp, beta)
            i_d, i_q = to_dq(*phase_currents(current_amp, beta))
            psi_d, psi_q = to_dq(*psi)
            if axis == "d":
                inductance = (psi_d - psi_d0) / i_d
            else:
                inductance = psi_q / i_q
            rows.append({
                "Current_A": current_amp, "Axis": axis, "Beta_deg": beta,
                "Id_A": i_d, "Iq_A": i_q,
                "PsiA_Wb": psi[0], "PsiB_Wb": psi[1], "PsiC_Wb": psi[2],
                "PsiD_Wb": psi_d, "PsiQ_Wb": psi_q,
                "L_H": inductance, "Torque_Nm": torque,
            })
            airgap.append((current_amp, axis, angles, br))
            print(f"I={current_amp} A {axis}-axis: psi_d = {psi_d:.4f} Wb, "
                  f"psi_q = {psi_q:.4f} Wb, L{axis} = {inductance * 1e3:.3f} mH, "
                  f"T = {torque:.2f} N*m")

    femm.closefemm()
    return rows, airgap


def inductances(rows):
    """current -> (Ld, Lq) in H."""
    table = {}
    for r in rows:
        table.setdefault(r["Current_A"], {})[r["Axis"]] = r["L_H"]
    return {i: (v["d"], v["q"]) for i, v in table.items()}


def save_results(rows, airgap):
    with open(RESULTS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with open(AIRGAP_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Current_A", "Axis", "Angle_deg", "Br_T"])
        for current_amp, axis, angles, br in airgap:
            writer.writerows((current_amp, axis, a, b) for a, b in zip(angles, br))


def load_results():
    """Rows and airgap B from a previous run's CSVs."""
    with open(RESULTS_FILE, newline="") as f:
        rows = [
            {k: (v if k == "Axis" else float(v)) for k, v in r.items()}
            for r in csv.DictReader(f)
        ]
    runs = {}
    with open(AIRGAP_FILE, newline="") as f:
        for r in csv.DictReader(f):
            angles, br = runs.setdefault((float(r["Current_A"]), r["Axis"]), ([], []))
            angles.append(float(r["Angle_deg"]))
            br.append(float(r["Br_T"]))
    airgap = [(i, axis, angles, br) for (i, axis), (angles, br) in runs.items()]
    return rows, airgap


def print_summary(rows, airgap):
    print(f"{'I [A]':>6} {'Ld [mH]':>8} {'Lq [mH]':>8} {'Lq/Ld':>6} "
          f"{'Ld-Lq [mH]':>11} {'ref Ld-Lq':>10}")
    for current_amp, (ld, lq) in inductances(rows).items():
        ref = REF_LD_MINUS_LQ.get(current_amp)
        ref_text = f"{ref * 1e3:>10.2f}" if ref is not None else ""
        print(f"{current_amp:>6.0f} {ld * 1e3:>8.3f} {lq * 1e3:>8.3f} {lq / ld:>6.2f} "
              f"{(ld - lq) * 1e3:>11.3f} {ref_text}")
    print(f"Reference unsaturated: Ld = {REF_LD * 1e3:.2f} mH, Lq = {REF_LQ * 1e3:.2f} mH")

    total = inductances(rows)
    print()
    print("Magnetizing (airgap flux) vs total, and the reference chart [mH]:")
    print(f"{'I [A]':>6} {'Lmd':>6} {'ref Ld':>7} {'diff':>6} {'Lmq':>6} {'ref Lq':>7} "
          f"{'diff':>6} {'Lsig_d':>7} {'Lsig_q':>7}")
    for current_amp, (lmd, lmq) in magnetizing_inductances(airgap).items():
        ld, lq = total[current_amp]
        ref_d = REF_CHART_LD.get(current_amp)
        ref_q = REF_CHART_LQ.get(current_amp)
        diff_d = f"{100 * (lmd - ref_d) / ref_d:>+5.0f}%" if ref_d else ""
        diff_q = f"{100 * (lmq - ref_q) / ref_q:>+5.0f}%" if ref_q else ""
        print(f"{current_amp:>6.0f} {lmd * 1e3:>6.2f} {ref_d * 1e3:>7.2f} {diff_d:>6} "
              f"{lmq * 1e3:>6.2f} {ref_q * 1e3:>7.2f} {diff_q:>6} "
              f"{(ld - lmd) * 1e3:>7.2f} {(lq - lmq) * 1e3:>7.2f}")
    print("Lsig = L - Lm: slot and tooth-tip leakage inductance.")


def plot_results(rows, airgap):
    table = inductances(rows)
    currents = list(table)
    ld = [table[i][0] * 1e3 for i in currents]
    lq = [table[i][1] * 1e3 for i in currents]

    fig, (ax_l, ax_diff) = plt.subplots(2, 1, sharex=True, figsize=(7, 7))
    ax_l.plot(currents, ld, "o-", label="Ld (FEMM)")
    ax_l.plot(currents, lq, "s-", label="Lq (FEMM)")
    ax_l.axhline(REF_LD * 1e3, color="C0", linestyle=":", label="Ld reference (unsaturated)")
    ax_l.axhline(REF_LQ * 1e3, color="C1", linestyle=":", label="Lq reference (unsaturated)")
    ax_l.set_ylabel("Inductance [mH]")
    ax_l.set_title("d- and q-axis Inductance vs Current"
                   + (" (magnets as air)" if MAGNETS_AS_AIR else ""))
    ax_l.grid(True)
    ax_l.legend(fontsize="small")

    ax_diff.plot(currents, [d - q for d, q in zip(ld, lq)], "o-", color="C2",
                 label="Ld - Lq (FEMM)")
    ax_diff.plot(list(REF_LD_MINUS_LQ), [v * 1e3 for v in REF_LD_MINUS_LQ.values()],
                 "x--", color="k", label="Ld - Lq (dq-model fit to test)")
    ax_diff.set_xlabel("Peak Current [A]")
    ax_diff.set_ylabel("Ld - Lq [mH]")
    ax_diff.grid(True)
    ax_diff.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(INDUCTANCE_PLOT_FILE)

    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(11, 4.5))
    for ax, axis in zip(axes, AXES):
        for current_amp, run_axis, angles, br in airgap:
            if run_axis == axis:
                ax.plot(angles, br, label=f"{current_amp} A")
        ax.set_title(f"Airgap radial B, current on the {axis}-axis")
        ax.set_xlabel("Mechanical angle [deg]")
        ax.grid(True)
    axes[0].set_ylabel("Br [T]")
    axes[1].legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(AIRGAP_PLOT_FILE)

    # Total vs magnetizing inductance vs the reference chart.
    magnetizing = magnetizing_inductances(airgap)
    fig, (ax_d, ax_q) = plt.subplots(1, 2, sharey=True, figsize=(11, 4.5))
    for ax, k, ref, name in ((ax_d, 0, REF_CHART_LD, "Ld"), (ax_q, 1, REF_CHART_LQ, "Lq")):
        ax.plot(currents, [table[i][k] * 1e3 for i in currents], "o-",
                label=f"{name} total (flux linkage)")
        ax.plot(currents, [magnetizing[i][k] * 1e3 for i in currents], "s-",
                label=f"{name} magnetizing (airgap flux)")
        ax.plot(list(ref), [v * 1e3 for v in ref.values()], "kx--", label=f"{name} reference")
        ax.set_title(f"{name}: total vs magnetizing vs reference")
        ax.set_xlabel("Peak Current [A]")
        ax.grid(True)
        ax.legend(fontsize="small")
    ax_d.set_ylabel("Inductance [mH]")
    fig.tight_layout()
    fig.savefig(MAGNETIZING_PLOT_FILE)

    plt.show()


if __name__ == "__main__":
    if "--plot" in sys.argv:
        rows, airgap = load_results()
    else:
        start_time = walltime.perf_counter()

        rows, airgap = run_inductance_sweep()
        save_results(rows, airgap)

        elapsed = walltime.perf_counter() - start_time
        print(f"Simulation run time: {elapsed:.1f} s")
    print_summary(rows, airgap)

    plot_results(rows, airgap)
