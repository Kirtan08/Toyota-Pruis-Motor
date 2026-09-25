import csv
import os
import shutil
import tempfile
import time as walltime

import cv2
import femm
import matplotlib.pyplot as plt

import config
import simulation

MODEL_FILE = "ToyotaPrius_MaxTorqueRotatingField.FEM"
RESULTS_FILE = "max_torque_rotating_field.csv"
TORQUE_PLOT_FILE = "max_torque_rotating_field.png"
TORQUE_ANGLE_PLOT_FILE = "max_torque_rotating_field_vs_angle.png"
CURRENT_PLOT_FILE = "max_torque_rotating_field_currents.png"
B_FIELD_VIDEO_FILE = "max_torque_rotating_field_b_field.mp4"
FLUX_LINES_VIDEO_FILE = "max_torque_rotating_field_flux_lines.mp4"


def _mo_showcontourplot(numcontours, al, au, ptype):
    """Workaround for a bug in the installed femm package: its
    mo_showcontourplot (site-packages/femm/__init__.py) builds the Lua
    call string with a stray comma instead of '+', which makes it pass 2
    arguments to callfemm() and raise a TypeError before ever reaching
    FEMM. Build the same command string by hand and call callfemm
    directly, bypassing the broken wrapper."""
    femm.callfemm(
        "mo_showcontourplot("
        + femm.numc(numcontours)
        + femm.numc(al)
        + femm.numc(au)
        + femm.quote(ptype)
        + ")"
    )


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

    frame_dir = tempfile.mkdtemp(prefix="max_torque_rotating_field_frames_")
    b_frames = []
    flux_frames = []

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

        # |B| density plot frame.
        femm.mo_showdensityplot(
            config.FieldPlotLegend,
            config.FieldPlotGrayscale,
            config.BPlotUpper,
            config.BPlotLower,
            "bmag",
        )
        b_frames.append(_capture_frame(frame_dir, step, "b"))

        # Flux-line (vector potential A contour) plot frame -- numcontours=-1
        # auto-ranges since the A range isn't known up front.
        _mo_showcontourplot(-1, 0, 0, "real")
        flux_frames.append(_capture_frame(frame_dir, step, "flux"))

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

    try:
        _frames_to_video(b_frames, B_FIELD_VIDEO_FILE, config.VideoFrameRate)
        _frames_to_video(flux_frames, FLUX_LINES_VIDEO_FILE, config.VideoFrameRate)
    finally:
        shutil.rmtree(frame_dir, ignore_errors=True)

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
