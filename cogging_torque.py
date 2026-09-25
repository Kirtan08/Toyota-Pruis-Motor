import csv
import os
import shutil
import tempfile
import time

import cv2
import femm
import matplotlib.pyplot as plt

import config
import simulation

OUTPUT_DIR = "Cogging_outputs"
BAND_BOUNDARY = simulation.SLIDING_BAND_NAME
MODEL_FILE = os.path.join(OUTPUT_DIR, "Cogging_outputs.FEM")
RESULTS_FILE = os.path.join(OUTPUT_DIR, "Cogging_outputs.csv")
PLOT_FILE = os.path.join(OUTPUT_DIR, "Cogging_outputs.png")
VIDEO_FILE = os.path.join(OUTPUT_DIR, "Cogging_outputs.mp4")
MESH_IMAGE_FILE = os.path.join(OUTPUT_DIR, "Cogging_mesh.bmp")


def _capture_frame(frame_dir, step, plot_type):
    """Zoom-to-fit the postprocessor view and save it as a numbered bitmap
    frame for later assembly into a video."""
    femm.mo_zoomnatural()
    frame_path = os.path.join(frame_dir, f"{plot_type}_{step:04d}.bmp")
    femm.mo_savebitmap(frame_path)
    return frame_path


def _frames_to_video(frame_paths, output_path, fps):
    """Assemble an ordered list of image frames into an .mp4 video."""
    first_frame = cv2.imread(frame_paths[0])
    height, width = first_frame.shape[:2]
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    writer = cv2.VideoWriter(output_path, fourcc, fps, (width, height))
    try:
        for frame_path in frame_paths:
            writer.write(cv2.imread(frame_path))
    finally:
        writer.release()


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

    frame_dir = tempfile.mkdtemp(prefix="cogging_frames_")
    b_frames = []

    total_steps = config.CoggingNumSteps + 1

    try:
        for step in range(total_steps):
            angle = step * config.CoggingStepAngle
            print(f"Step {step + 1}/{total_steps}: angle = {angle:.4f} deg")

            femm.mi_modifyboundprop(BAND_BOUNDARY, 10, angle)

            femm.mi_analyze(1)
            femm.mi_loadsolution()

            torque = femm.mo_gapintegral(BAND_BOUNDARY, 0)
            angles.append(angle)
            torques.append(torque)

            # |B| density plot frame.
            femm.mo_showdensityplot(
                config.FieldPlotLegend,
                config.FieldPlotGrayscale,
                config.BPlotUpper,
                config.BPlotLower,
                "bmag",
            )
            b_frames.append(_capture_frame(frame_dir, step, "b"))

            femm.mo_close()

        _frames_to_video(b_frames, VIDEO_FILE, config.VideoFrameRate)
    finally:
        shutil.rmtree(frame_dir, ignore_errors=True)

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
