import csv
import os
import sys
import time as walltime

import cv2
import femm
import matplotlib.pyplot as plt

import config
import simulation

if len(sys.argv) > 1:
    config.TorqueVsCurrentAmps = [float(sys.argv[1])]

RESULTS_FILE = "torque_vs_current.csv"
KT_RESULTS_FILE = "torque_vs_current_kt.csv"
TORQUE_PHASE_PLOT_FILE = "torque_vs_current_phase.png"
KT_PLOT_FILE = "torque_vs_current_kt.png"
MODEL_FILE = "ToyotaPrius_TorqueVsCurrent.FEM"
FRAMES_DIR = "Locked_rotor_frames"
VIDEO_FILE = "Locked_rotor_B.mp4"
VIDEO_FPS = 4

# |B| density plot scale, T -- fixed so every frame shares the same colors.
B_PLOT_MIN = 0
B_PLOT_MAX = 2


def save_flux_density_frame(path):
    """Save the loaded solution's |B| density plot, zoomed to fit."""
    femm.mo_zoomnatural()
    femm.mo_hidepoints()
    femm.mo_showdensityplot(1, 0, B_PLOT_MAX, B_PLOT_MIN, "bmag")
    femm.mo_savebitmap(path)


def run_torque_vs_current_sweep():
    os.makedirs(FRAMES_DIR, exist_ok=True)

    femm.openfemm()
    simulation.build_model()

    simulation.pause("Geometry, materials, and boundary conditions complete.")


    femm.mi_modifyboundprop(
        simulation.SLIDING_BAND_NAME, 10,
        config.MaxTorqueInitialAngle,
    )

    all_currents = []
    all_phases = []
    all_torques = []
    peak_currents = []
    peak_phases = []
    peak_torques = []
    frames = []

    for current_amp in config.TorqueVsCurrentAmps:
        phases = []
        torques = []

        for step in range(config.TorqueVsCurrentPhaseSteps + 1):
            phase = config.TorqueVsCurrentPhaseInit + step * config.TorqueVsCurrentPhaseStep

            i_a = current_amp * simulation.sind(phase)
            i_b = current_amp * simulation.sind(phase + 120)
            i_c = current_amp * simulation.sind(phase + 240)
            femm.mi_modifycircprop("A", 1, i_a)
            femm.mi_modifycircprop("B", 1, i_b)
            femm.mi_modifycircprop("C", 1, i_c)

            femm.mi_saveas(MODEL_FILE)
            femm.mi_createmesh()
            femm.mi_analyze(0)
            femm.mi_loadsolution()

            torque = femm.mo_gapintegral(simulation.SLIDING_BAND_NAME, 0)

            frame_file = os.path.join(FRAMES_DIR, f"B_{current_amp:03.0f}A_{step:02d}.bmp")
            save_flux_density_frame(frame_file)
            frames.append((frame_file, current_amp, phase, torque))
            femm.mo_close()

            phases.append(phase)
            torques.append(torque)

            print(f"I={current_amp} A :: phase step {step} :: {config.TorqueVsCurrentPhaseSteps}")

        peak_torque = max(torques)
        peak_phase = phases[torques.index(peak_torque)]

        all_currents.extend([current_amp] * len(phases))
        all_phases.extend(phases)
        all_torques.extend(torques)
        peak_currents.append(current_amp)
        peak_phases.append(peak_phase)
        peak_torques.append(peak_torque)

    femm.closefemm()
    return {
        "currents": all_currents, "phases": all_phases, "torques": all_torques,
        "peak_currents": peak_currents, "peak_phases": peak_phases, "peak_torques": peak_torques,
        "frames": frames,
    }


def make_video(frames):
    """Stitch the saved |B| frames into an MP4, labeling each with its
    current, current angle and torque."""
    if not frames:
        return
    height, width = cv2.imread(frames[0][0]).shape[:2]
    writer = cv2.VideoWriter(
        VIDEO_FILE, cv2.VideoWriter_fourcc(*"mp4v"), VIDEO_FPS, (width, height)
    )
    for path, current_amp, phase, torque in frames:
        img = cv2.imread(path)
        if img.shape[:2] != (height, width):
            img = cv2.resize(img, (width, height))
        label = f"I = {current_amp:.0f} A   angle = {phase:.0f} deg   T = {torque:.1f} N*m"
        cv2.putText(img, label, (10, height - 15), cv2.FONT_HERSHEY_SIMPLEX,
                    0.6, (0, 0, 0), 2, cv2.LINE_AA)
        writer.write(img)
    writer.release()


def save_results(results):
    with open(RESULTS_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Current_A", "Phase_deg", "Torque_Nm"])
        writer.writerows(zip(results["currents"], results["phases"], results["torques"]))

    with open(KT_RESULTS_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Current_A", "PeakPhase_deg", "PeakTorque_Nm"])
        writer.writerows(zip(results["peak_currents"], results["peak_phases"], results["peak_torques"]))


def plot_results(results):
    plt.figure()
    for current_amp in config.TorqueVsCurrentAmps:
        phases = [p for c, p in zip(results["currents"], results["phases"]) if c == current_amp]
        torques = [t for c, t in zip(results["currents"], results["torques"]) if c == current_amp]
        plt.plot(phases, torques, ".-", label=f"{current_amp} A")
    plt.xlabel("Electrical Angle, deg")
    plt.ylabel("Torque, N*m")
    plt.legend()
    plt.grid(True)
    plt.title("Torque vs Current-Phase, per Current Amplitude")
    plt.savefig(TORQUE_PHASE_PLOT_FILE)

    plt.figure()
    plt.plot(results["peak_currents"], results["peak_torques"], ".-", color="g")
    plt.xlabel("Peak Current, A")
    plt.ylabel("Peak Torque, N*m")
    plt.grid(True)
    plt.title("Peak Torque vs Peak Current (Kt curve)")
    plt.savefig(KT_PLOT_FILE)

    plt.show()


if __name__ == "__main__":
    start_time = walltime.perf_counter()

    results = run_torque_vs_current_sweep()
    save_results(results)
    make_video(results["frames"])

    elapsed = walltime.perf_counter() - start_time
    print(f"Simulation run time: {elapsed:.1f} s")
    print(f"Peak torques by current: {list(zip(results['peak_currents'], results['peak_torques']))}")

    plot_results(results)
