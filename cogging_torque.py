import csv
import os
import time

import femm
import matplotlib.pyplot as plt

import config
import simulation

OUTPUT_DIR = "Cogging_outputs"
BAND_BOUNDARY = simulation.SLIDING_BAND_NAME
MODEL_FILE = os.path.join(OUTPUT_DIR, "Cogging_outputs.FEM")
RESULTS_FILE = os.path.join(OUTPUT_DIR, "Cogging_outputs.csv")
PLOT_FILE = os.path.join(OUTPUT_DIR, "Cogging_outputs.png")
MESH_IMAGE_FILE = os.path.join(OUTPUT_DIR, "Cogging_mesh.bmp")


def run_cogging_sweep():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    femm.openfemm()
    # Cogging torque is measured unexcited, so there's no winding to carry
    # current: windings=False fills the slots with plain (coarsely meshed)
    # Air instead. split_stator_mesh/split_rotor_mesh give the tooth tips
    # and the rotor edge facing the stator (both torque-critical, like the
    # airgap) a dense mesh while the slot bodies and bulk yoke/rotor steel
    # stay coarse.
    simulation.build_model(windings=False, split_stator_mesh=True, split_rotor_mesh=True)

    # The one pause in the whole run: inspect the geometry, materials, and
    # boundary conditions before committing to the (potentially long) sweep.
    simulation.pause("Geometry, materials, and boundary conditions complete.")

    femm.mi_saveas(MODEL_FILE)

    # Automatic mesh generation, as in reference.txt -- smartmesh (enabled in
    # simulation.build_model) sizes the mesh from the geometry, this just
    # triggers it once up front rather than implicitly on the first analyze.
    femm.mi_createmesh()

    # Snapshot of the initial (step 0) mesh only, zoomed to fit, so the
    # mesh density can be inspected without re-opening the model.
    femm.mi_showmesh()
    femm.mi_zoomnatural()
    femm.mi_savebitmap(MESH_IMAGE_FILE)

    angles = []
    torques = []

    total_steps = config.CoggingNumSteps + 1

    for step in range(total_steps):
        angle = step * config.CoggingStepAngle
        print(f"Step {step + 1}/{total_steps}: angle = {angle:.4f} deg")

        femm.mi_modifyboundprop(BAND_BOUNDARY, 10, angle)

        femm.mi_analyze(1)
        femm.mi_loadsolution()

        torque = femm.mo_gapintegral(BAND_BOUNDARY, 0)
        angles.append(angle)
        torques.append(torque)

        femm.mo_close()

    femm.closefemm()
    return angles, torques


def save_results(angles, torques, path=RESULTS_FILE):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["MechanicalAngle_deg", "Torque_Nm"])
        writer.writerows(zip(angles, torques))


def plot_results(angles, torques, path=PLOT_FILE):
    plt.figure()
    plt.plot(angles, torques, marker="o")
    plt.xlabel("Rotor Mechanical Angle, deg")
    plt.ylabel("Cogging Torque, N*m")
    plt.grid(True)
    plt.title(
        f"Cogging Torque -- {config.Nslots} slots / {config.Npoles} poles "
        f"(period = {config.CoggingPeriodAngle:g} deg)"
    )
    plt.savefig(path)
    plt.show()


if __name__ == "__main__":
    start_time = time.perf_counter()

    angles, torques = run_cogging_sweep()
    save_results(angles, torques)

    elapsed = time.perf_counter() - start_time
    print(f"Simulation run time: {elapsed:.1f} s")

    plot_results(angles, torques)
