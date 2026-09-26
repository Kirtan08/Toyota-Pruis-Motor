import csv
import time as walltime

import femm
import matplotlib.pyplot as plt

import config
import simulation

MODEL_FILE = "ToyotaPrius_StaticRotorRotatingCurrent.FEM"
RESULTS_FILE = "static_rotor_rotating_current.csv"
TORQUE_TIME_PLOT_FILE = "static_rotor_rotating_current_vs_time.png"

# Current waveform
FREQ = config.Freq          # Hz, 66.7 Hz (8 poles at 1000 RPM)
SIM_TIME = 15e-3            # s, total simulated time (one electrical period)
STEP_TIME = 1e-3            # s, time step
NUM_STEPS = int(round(SIM_TIME / STEP_TIME))
PHASE = 0                   # deg, phase A current angle at t=0


def phase_currents(t):
    """3-phase currents at time t (s): config.Current amplitude, positive
    (A-B-C) sequence with B and C lagging A by 120 and 240 deg. At t=0 this
    gives IA = 0, IB = -8.66 A, IC = 8.66 A (for 10 A peak)."""
    return tuple(
        config.Current * simulation.sind(360 * FREQ * t + PHASE - shift)
        for shift in (0, 120, 240)
    )


def run_static_rotor_sweep():
    """Static-rotor sweep: hold the rotor fixed at MaxTorqueInitialAngle and
    advance the 3-phase currents in time from 0 to SIM_TIME in STEP_TIME
    steps, tracing torque vs time for a rotating stator field."""
    femm.openfemm()
    simulation.build_model()

    simulation.pause("Geometry, materials, and boundary conditions complete.")

    # Rotor position is set once and never touched again.
    femm.mi_modifyboundprop(
        simulation.SLIDING_BAND_NAME, 10,
        config.MaxTorqueInitialAngle - config.SectorAngle,
    )

    times = []
    torques = []
    currents = []

    for step in range(NUM_STEPS + 1):
        t = step * STEP_TIME

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
        torques.append(torque)
        currents.append((i_a, i_b, i_c))

        print(
            f"{step} :: {NUM_STEPS}  t = {t * 1e3:.1f} ms  "
            f"IA = {i_a:.3f} A, IB = {i_b:.3f} A, IC = {i_c:.3f} A"
        )

    femm.closefemm()
    return times, torques, currents


def save_results(times, torques, currents, path=RESULTS_FILE):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Time_ms", "Torque_Nm", "IA_A", "IB_A", "IC_A"])
        writer.writerows(
            (t * 1e3, torque, i_a, i_b, i_c)
            for t, torque, (i_a, i_b, i_c) in zip(times, torques, currents)
        )


def plot_results(times, torques, currents):
    times_ms = [t * 1e3 for t in times]

    fig, (ax_torque, ax_current) = plt.subplots(2, 1, sharex=True)

    ax_torque.plot(times_ms, torques, ".-", color="b")
    ax_torque.set_ylabel("Torque, N*m")
    ax_torque.grid(True)
    ax_torque.set_title(
        f"Static Rotor, Rotating Current -- Torque vs Time\n"
        f"({config.Current:g} A peak, {FREQ:.1f} Hz)"
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

    times, torques, currents = run_static_rotor_sweep()
    save_results(times, torques, currents)

    elapsed = walltime.perf_counter() - start_time
    print(f"Simulation run time: {elapsed:.1f} s")
    print(f"Max torque: {max(torques):.2f} N*m at {times[torques.index(max(torques))] * 1e3:.1f} ms")

    plot_results(times, torques, currents)
