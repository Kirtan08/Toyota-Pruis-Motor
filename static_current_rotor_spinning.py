import csv
import time as walltime

import femm
import matplotlib.pyplot as plt

import config
import simulation

MODEL_FILE = "ToyotaPrius_StaticCurrentRotorSpinning.FEM"
RESULTS_FILE = "static_current_rotor_spinning.csv"
TORQUE_ANGLE_PLOT_FILE = "static_current_rotor_spinning_vs_angle.png"

# Fixed phase currents, A
I_A = 0
I_B = -8.660
I_C = 8.660


def run_static_current_sweep():
    """Static-current sweep: hold the 3-phase currents fixed at I_A/I_B/I_C
    and rotate ONLY the rotor (via the sliding band) across
    MaxTorqueStaticSweepAngle (90 mech deg = one electrical period), tracing
    the torque-angle curve for a single fixed current vector."""
    femm.openfemm()
    simulation.build_model()

    simulation.pause("Geometry, materials, and boundary conditions complete.")

    # Set once and never touched again.
    i_a, i_b, i_c = I_A, I_B, I_C
    femm.mi_modifycircprop("A", 1, i_a)
    femm.mi_modifycircprop("B", 1, i_b)
    femm.mi_modifycircprop("C", 1, i_c)
    print(f"Fixed currents: IA = {i_a:.3f} A, IB = {i_b:.3f} A, IC = {i_c:.3f} A")

    angles = []
    torques = []

    for step in range(config.MaxTorqueStaticNumSteps + 1):
        rotation = step * config.MaxTorqueStaticStepAngle
        angle = config.MaxTorqueInitialAngle + rotation

        # Same sector-vs-full-model rotor reference offset as
        # max_torque.py's run_max_torque_sweep().
        femm.mi_modifyboundprop(
            simulation.SLIDING_BAND_NAME, 10, angle - config.SectorAngle
        )

        femm.mi_saveas(MODEL_FILE)
        femm.mi_createmesh()
        femm.mi_analyze(0)
        femm.mi_loadsolution()

        torque = femm.mo_gapintegral(simulation.SLIDING_BAND_NAME, 0)
        femm.mo_close()

        angles.append(rotation)
        torques.append(torque)

        print(f"{step} :: {config.MaxTorqueStaticNumSteps}  rotation = {rotation:.2f} deg")

    femm.closefemm()
    return angles, torques, (i_a, i_b, i_c)


def save_results(angles, torques, currents, path=RESULTS_FILE):
    i_a, i_b, i_c = currents
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["RotorRotation_deg", "Torque_Nm", "IA_A", "IB_A", "IC_A"])
        writer.writerows((a, t, i_a, i_b, i_c) for a, t in zip(angles, torques))


def plot_results(angles, torques, currents):
    i_a, i_b, i_c = currents
    plt.figure()
    plt.plot(angles, torques, ".-", color="b")
    plt.xlabel("Rotor Rotation, mech deg")
    plt.ylabel("Torque, N*m")
    plt.grid(True)
    plt.title(
        f"Static Current, Rotor Spinning -- Torque vs Rotor Position\n"
        f"(IA = {i_a:.2f} A, IB = {i_b:.2f} A, IC = {i_c:.2f} A)"
    )
    plt.savefig(TORQUE_ANGLE_PLOT_FILE)

    plt.show()


if __name__ == "__main__":
    start_time = walltime.perf_counter()

    angles, torques, currents = run_static_current_sweep()
    save_results(angles, torques, currents)

    elapsed = walltime.perf_counter() - start_time
    print(f"Simulation run time: {elapsed:.1f} s")
    print(f"Max torque: {max(torques):.2f} N*m at {angles[torques.index(max(torques))]:.2f} deg")

    plot_results(angles, torques, currents)
