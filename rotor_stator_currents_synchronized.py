import csv
import time as walltime

import femm
import matplotlib.pyplot as plt

import config
import simulation

MODEL_FILE = "ToyotaPrius_RotorStatorCurrentsSynchronized.FEM"
RESULTS_FILE = "rotor_stator_currents_synchronized.csv"
TORQUE_TIME_PLOT_FILE = "rotor_stator_currents_synchronized_vs_time.png"

# Current waveform
FREQ = config.Freq          # Hz, 66.7 Hz (8 poles at 1000 RPM)
SIM_TIME = 15e-3            # s, total simulated time (one electrical period)
STEP_TIME = 1e-3            # s, time step
NUM_STEPS = int(round(SIM_TIME / STEP_TIME))
PHASE = config.Phase        # deg, phase A current angle at t=0 (120 deg)

# Rotor motion
SPEED_RPM = 1000            # RPM, synchronous speed for FREQ and config.Npoles
ROTOR_SPEED = SPEED_RPM * 360 / 60   # mech deg/s -- 6 deg per 1 ms step
# Rotor turns in the POSITIVE sliding-band angle direction, matching the
# reference script's positive mi_moverotate steps. This keeps it synchronized
# with the A-C-B field: a -1 run gave torque swinging between about +/-219
# N*m instead of staying constant.
ROTOR_DIRECTION = 1


def phase_currents(t):
    """3-phase currents at time t (s), the reference script's waveform:
    config.Current amplitude, A-C-B sequence with B and C shifted +120 and
    +240 deg from A. At t=0 this gives IA = 8.66 A, IB = -8.66 A, IC = 0
    (for 10 A peak, 120 deg phase)."""
    return tuple(
        config.Current * simulation.sind(360 * FREQ * t + PHASE + shift)
        for shift in (0, 120, 240)
    )


def rotor_rotation(t):
    """Rotor rotation (mech deg) from its initial position at time t (s)."""
    return ROTOR_DIRECTION * ROTOR_SPEED * t


def run_synchronized_sweep():
    """Synchronized sweep: advance the 3-phase currents in time from 0 to
    SIM_TIME in STEP_TIME steps while rotating the rotor (via the sliding
    band) at SPEED_RPM from MaxTorqueInitialAngle, tracing torque vs time
    with the rotor turning together with the stator field."""
    femm.openfemm()
    simulation.build_model()

    simulation.pause("Geometry, materials, and boundary conditions complete.")

    times = []
    rotations = []
    torques = []
    currents = []

    for step in range(NUM_STEPS + 1):
        t = step * STEP_TIME

        # No "- config.SectorAngle" offset here, unlike the other scripts:
        # the sector's magnets are drawn in the same position as the
        # reference's full model, so the rotor starts at the reference's own
        # +7.5 deg. With anti-periodic boundaries, the extra -45 deg (one
        # pole) reversed the magnet polarity: -56 N*m at t=0 instead of
        # about +210 N*m.
        rotation = rotor_rotation(t)
        angle = config.MaxTorqueInitialAngle + rotation
        femm.mi_modifyboundprop(simulation.SLIDING_BAND_NAME, 10, angle)

        i_a, i_b, i_c = phase_currents(t)
        femm.mi_modifycircprop("A", 1, i_a)
        femm.mi_modifycircprop("B", 1, i_b)
        femm.mi_modifycircprop("C", 1, i_c)

        femm.mi_saveas(MODEL_FILE)
        femm.mi_createmesh()
        femm.mi_analyze(0)
        femm.mi_loadsolution()

        torque = femm.mo_gapintegral(simulation.SLIDING_BAND_NAME, 0)
        femm.mo_close()

        times.append(t)
        rotations.append(rotation)
        torques.append(torque)
        currents.append((i_a, i_b, i_c))

        print(
            f"{step} :: {NUM_STEPS}  t = {t * 1e3:.1f} ms  rotation = {rotation:.2f} deg  "
            f"IA = {i_a:.3f} A, IB = {i_b:.3f} A, IC = {i_c:.3f} A"
        )

    femm.closefemm()
    return times, rotations, torques, currents


def save_results(times, rotations, torques, currents, path=RESULTS_FILE):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Time_ms", "RotorRotation_deg", "Torque_Nm", "IA_A", "IB_A", "IC_A"])
        writer.writerows(
            (t * 1e3, rotation, torque, i_a, i_b, i_c)
            for t, rotation, torque, (i_a, i_b, i_c) in zip(times, rotations, torques, currents)
        )


def plot_results(times, torques, currents):
    times_ms = [t * 1e3 for t in times]

    fig, (ax_torque, ax_current) = plt.subplots(2, 1, sharex=True)

    ax_torque.plot(times_ms, torques, ".-", color="b")
    ax_torque.set_ylabel("Torque, N*m")
    ax_torque.grid(True)
    ax_torque.set_title(
        f"Rotor and Stator Currents Synchronized -- Torque vs Time\n"
        f"({config.Current:g} A peak, {FREQ:.1f} Hz, rotor {SPEED_RPM} RPM)"
    )

    for k, name in enumerate(("IA", "IB", "IC")):
        ax_current.plot(times_ms, [c[k] for c in currents], ".-", label=name)
    ax_current.set_xlabel("Time, ms")
    ax_current.set_ylabel("Current, A")
    ax_current.legend()
    ax_current.grid(True)

    fig.savefig(TORQUE_TIME_PLOT_FILE)

    plt.show()


if __name__ == "__main__":
    start_time = walltime.perf_counter()

    times, rotations, torques, currents = run_synchronized_sweep()
    save_results(times, rotations, torques, currents)

    elapsed = walltime.perf_counter() - start_time
    print(f"Simulation run time: {elapsed:.1f} s")
    print(f"Max torque: {max(torques):.2f} N*m at {times[torques.index(max(torques))] * 1e3:.1f} ms")
    print(f"Mean torque: {sum(torques) / len(torques):.2f} N*m")

    plot_results(times, torques, currents)
