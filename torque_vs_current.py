import csv
import sys
import time as walltime

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


def run_torque_vs_current_sweep():

    femm.openfemm()
    simulation.build_model()

    simulation.pause("Geometry, materials, and boundary conditions complete.")

    # Rotor position is set once and never touched again.
    femm.mi_modifyboundprop(
        simulation.SLIDING_BAND_NAME, 10,
        config.MaxTorqueInitialAngle - config.SectorAngle + 30,
    )

    all_currents = []
    all_phases = []
    all_torques = []
    peak_currents = []
    peak_phases = []
    peak_torques = []

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
    }


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

    elapsed = walltime.perf_counter() - start_time
    print(f"Simulation run time: {elapsed:.1f} s")
    print(f"Peak torques by current: {list(zip(results['peak_currents'], results['peak_torques']))}")

    plot_results(results)
