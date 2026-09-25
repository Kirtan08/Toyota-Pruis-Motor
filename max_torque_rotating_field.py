import csv
import time as walltime

import femm
import matplotlib.pyplot as plt

import config
import simulation

MODEL_FILE = "ToyotaPrius_MaxTorqueRotatingField.FEM"
RESULTS_FILE = "max_torque_rotating_field.csv"
TORQUE_PLOT_FILE = "max_torque_rotating_field.png"
TORQUE_ANGLE_PLOT_FILE = "max_torque_rotating_field_vs_angle.png"
CURRENT_PLOT_FILE = "max_torque_rotating_field_currents.png"


def run_max_torque_rotating_field_sweep():
    """Synchronized sweep, per prius_motor_full_model.py's actual method:
    the rotor (via the sliding band) and the 3-phase currents advance
    together in time, over one electrical period -- same calculation as
    max_torque.py's run_max_torque_sweep(), just at the coarser 15-step
    resolution the reference blog used. Contrast with max_torque_static.py
    (currents fixed, rotor rotates)."""
    femm.openfemm()
    simulation.build_model()

    simulation.pause("Geometry, materials, and boundary conditions complete.")

    angle = config.MaxTorqueInitialAngle
    sim_time = 0.0

    times = []
    angles = []
    torques = []
    currents_a, currents_b, currents_c = [], [], []

    for step in range(config.MaxTorqueRotatingFieldNumSteps + 1):
        # Same sector-vs-full-model rotor reference offset as
        # max_torque.py's run_max_torque_sweep().
        femm.mi_modifyboundprop(
            simulation.SLIDING_BAND_NAME, 10, angle - config.SectorAngle
        )

        i_a = config.Current * simulation.sind(360 * config.Freq * sim_time + config.Phase)
        i_b = config.Current * simulation.sind(360 * config.Freq * sim_time + config.Phase + 120)
        i_c = config.Current * simulation.sind(360 * config.Freq * sim_time + config.Phase + 240)
        femm.mi_modifycircprop("A", 1, i_a)
        femm.mi_modifycircprop("B", 1, i_b)
        femm.mi_modifycircprop("C", 1, i_c)

        femm.mi_saveas(MODEL_FILE)
        femm.mi_createmesh()
        femm.mi_analyze(0)
        femm.mi_loadsolution()

        torque = femm.mo_gapintegral(simulation.SLIDING_BAND_NAME, 0)
        femm.mo_close()

        times.append(sim_time)
        angles.append(angle)
        torques.append(torque)
        currents_a.append(i_a)
        currents_b.append(i_b)
        currents_c.append(i_c)

        print(f"{step} :: {config.MaxTorqueRotatingFieldNumSteps}")

        sim_time += config.MaxTorqueRotatingFieldStepTime
        angle += config.MaxTorqueRotatingFieldStepAngle

    femm.closefemm()
    return times, angles, torques, currents_a, currents_b, currents_c


def save_results(times, angles, torques, currents_a, currents_b, currents_c, path=RESULTS_FILE):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["Time_s", "MechanicalAngle_deg", "Torque_Nm", "IA_A", "IB_A", "IC_A"]
        )
        writer.writerows(zip(times, angles, torques, currents_a, currents_b, currents_c))


def plot_results(times, angles, torques, currents_a, currents_b, currents_c):
    plt.figure()
    plt.plot(times, currents_a, ".-", label="Phase A")
    plt.plot(times, currents_b, ".-", label="Phase B")
    plt.plot(times, currents_c, ".-", label="Phase C")
    plt.xlabel("Time, s")
    plt.ylabel("Current, A")
    plt.legend()
    plt.grid(True)
    plt.title("Synchronized Sweep (15 steps) -- Phase Currents")
    plt.savefig(CURRENT_PLOT_FILE)

    plt.figure()
    plt.plot(times, torques, ".-", color="b")
    plt.xlabel("Time, s")
    plt.ylabel("Torque, N*m")
    plt.grid(True)
    plt.title(
        f"Synchronized Sweep (15 steps) -- Torque ({config.Current} A, "
        f"{config.SpeedRPM} RPM, {config.Freq:.2f} Hz)"
    )
    plt.savefig(TORQUE_PLOT_FILE)

    plt.figure()
    plt.plot(angles, torques, ".-", color="b")
    plt.xlabel("Rotor Mechanical Angle, deg")
    plt.ylabel("Torque, N*m")
    plt.grid(True)
    plt.title(
        f"Synchronized Sweep (15 steps) -- Torque vs Rotor Position "
        f"({config.Current} A, {config.SpeedRPM} RPM, {config.Freq:.2f} Hz)"
    )
    plt.savefig(TORQUE_ANGLE_PLOT_FILE)

    plt.show()


if __name__ == "__main__":
    start_time = walltime.perf_counter()

    times, angles, torques, ia, ib, ic = run_max_torque_rotating_field_sweep()
    save_results(times, angles, torques, ia, ib, ic)

    elapsed = walltime.perf_counter() - start_time
    print(f"Simulation run time: {elapsed:.1f} s")
    print(f"Max torque: {max(torques):.2f} N*m")

    plot_results(times, angles, torques, ia, ib, ic)
