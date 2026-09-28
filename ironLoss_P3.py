"""Iron loss (part 3): FFT of the stator B waveforms and iron loss map.

Takes the tooth and back-iron B waveforms from ironLoss_P2.py (one
electrical cycle at 1000 RPM, per current at its MTPA angle), splits them
into harmonics with an FFT, and applies the loss-separation model fitted in
ironLoss_P1.py to each odd harmonic q = 1, 3, 5, 7, 9:
    P = sum_q  Kh*(q*f0)*Bq^2 + Kc*(q*f0*Bq)^2 + Ke*(q*f0*Bq)^1.5   [W/kg]
The B waveform shape is assumed not to change with speed (magnetostatic
solve, MTPA current angle at every speed), so only f0 scales with RPM.

Two ways to get the loss in W:
  - reference: tooth B applied to the whole stator core mass, as in the
    reference study (which notes this is not correct);
  - split: tooth B over the tooth mass, back-iron B over the back-iron mass.

FFT: the P2 sweep covers exactly one electrical period, so after dropping
its last sample (a repeat of the first) FFT bin k is exactly harmonic k --
no need to tile the waveform as the reference did.

Usage:
    python ironLoss_P3.py
"""
import math

import matplotlib.pyplot as plt
import numpy as np

import config
from ironLoss_P1 import KC, KE, KH, SF
from ironLoss_P2 import POST_HEIGHT, by_current, load_results
from magnet_operating_point import FREQ

PLOT_PREFIX = "ironLoss_P3_"
RESULTS_FILE = "ironLoss_P3.csv"

HARMONICS = [1, 3, 5, 7, 9]
RPMS = [200, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000]

# Loss coefficients per kg of steel, corrected for the stacking factor as in
# the reference (Kh/0.97, Kc/0.97, Ke/sqrt(0.97)).
KH_SF = KH / SF
KC_SF = KC / SF
KE_SF = KE / math.sqrt(SF)

DENSITY = 7650               # kg/m^3
CORE_VOLUME = 0.00244055     # m^3, whole stator core (teeth + back iron), reference value
IN = 0.0254                  # m per in
# Back iron: annulus from the slot bottoms to the stator OD over the stack.
BACKIRON_VOLUME = (math.pi * ((config.StatorOD / 2) ** 2
                              - (config.StatorID / 2 + POST_HEIGHT) ** 2)
                   * IN ** 2 * config.Length * IN)
TEETH_VOLUME = CORE_VOLUME - BACKIRON_VOLUME

# Reference tooth B harmonic amplitudes [T] (rows: currents, columns: q).
REFERENCE_B_TOOTH = {
    0: [0.94031, 0.02239, 0.14721, 0.06210, 0.00912],
    50: [1.43650, 0.07363, 0.14686, 0.04078, 0.01904],
    75: [1.60251, 0.16855, 0.11853, 0.02913, 0.00592],
    100: [1.70893, 0.24756, 0.11589, 0.02463, 0.00889],
    125: [1.77942, 0.30065, 0.13350, 0.02036, 0.00886],
    150: [1.82498, 0.33113, 0.15539, 0.02707, 0.01328],
    200: [1.90894, 0.35532, 0.15714, 0.04976, 0.01947],
    250: [1.96453, 0.36682, 0.15219, 0.06258, 0.01975],
}


def spectrum(b):
    """Amplitude [T] and phase [rad] of every harmonic k = 0..N/2 of one
    period of samples b, so that b(t) = sum_k amp_k*cos(2*pi*k*f0*t + phase_k)."""
    fts = np.fft.rfft(b) / len(b)
    amp = np.abs(fts)
    amp[1:] *= 2
    if len(b) % 2 == 0:
        amp[-1] /= 2   # Nyquist bin is not doubled
    return amp, np.angle(fts)


def harmonics(rows):
    """{current: {'tooth'|'backiron': (times_s, b, amp, phase)}} over one
    period (the repeated last sample dropped)."""
    out = {}
    for current_amp, group in by_current(rows).items():
        group = group[:-1]
        times = np.array([r["Time_ms"] for r in group]) / 1e3
        out[current_amp] = {}
        for region, key in (("tooth", "Btooth_T"), ("backiron", "Bbackiron_T")):
            b = np.array([r[key] for r in group])
            out[current_amp][region] = (times, b, *spectrum(b))
    return out


def specific_loss(amp, rpm):
    """Iron loss [W/kg] at `rpm` from harmonic amplitudes amp[q]."""
    f0 = config.Npoles * rpm / 120
    loss = 0.0
    for q in HARMONICS:
        fb = q * f0 * amp[q]
        loss += KH_SF * q * f0 * amp[q] ** 2 + KC_SF * fb ** 2 + KE_SF * fb ** 1.5
    return loss


def loss_maps(spectra):
    """Iron loss [W] per (current, rpm): reference (tooth B over the whole
    core) and split (tooth B over the teeth, back-iron B over the back iron)."""
    currents = list(spectra)
    ref = np.zeros((len(currents), len(RPMS)))
    split = np.zeros_like(ref)
    for i, current_amp in enumerate(currents):
        tooth = spectra[current_amp]["tooth"][2]
        back = spectra[current_amp]["backiron"][2]
        for j, rpm in enumerate(RPMS):
            p_tooth = specific_loss(tooth, rpm)
            p_back = specific_loss(back, rpm)
            ref[i, j] = p_tooth * CORE_VOLUME * DENSITY
            split[i, j] = (p_tooth * TEETH_VOLUME + p_back * BACKIRON_VOLUME) * DENSITY
    return currents, ref, split


def print_harmonics(spectra):
    print(f"Fundamental f0 = {FREQ:.1f} Hz at 1000 RPM")
    print(f"Volumes: core {CORE_VOLUME * 1e6:.0f} cm^3, teeth {TEETH_VOLUME * 1e6:.0f} cm^3, "
          f"back iron {BACKIRON_VOLUME * 1e6:.0f} cm^3")
    head = " ".join(f"{'B' + str(q):>7}" for q in HARMONICS)
    for region in ("tooth", "backiron"):
        print(f"\n{region} B harmonic amplitudes [T]"
              + ("  (reference B1 in brackets)" if region == "tooth" else ""))
        print(f"{'I [A]':>6} {head}")
        for current_amp, s in spectra.items():
            amp = s[region][2]
            line = f"{current_amp:>6.0f} " + " ".join(f"{amp[q]:>7.4f}" for q in HARMONICS)
            if region == "tooth" and current_amp in REFERENCE_B_TOOTH:
                line += f"   ({REFERENCE_B_TOOTH[current_amp][0]:.4f})"
            print(line)


def print_losses(currents, ref, split):
    for name, table in (("reference: tooth B x whole core", ref),
                        ("split: teeth + back iron", split)):
        print(f"\nIron loss [W], {name}")
        print(f"{'I [A]':>6} " + " ".join(f"{rpm:>6}" for rpm in RPMS))
        for current_amp, row in zip(currents, table):
            print(f"{current_amp:>6.0f} " + " ".join(f"{p:>6.0f}" for p in row))


def save_results(spectra, currents, ref, split):
    with open(RESULTS_FILE, "w") as f:
        f.write("Current_A,RPM,"
                + ",".join(f"Btooth{q}_T" for q in HARMONICS) + ","
                + ",".join(f"Bbackiron{q}_T" for q in HARMONICS)
                + ",PironRef_W,PironSplit_W\n")
        for i, current_amp in enumerate(currents):
            tooth = spectra[current_amp]["tooth"][2]
            back = spectra[current_amp]["backiron"][2]
            for j, rpm in enumerate(RPMS):
                f.write(f"{current_amp:g},{rpm},"
                        + ",".join(f"{tooth[q]:.5f}" for q in HARMONICS) + ","
                        + ",".join(f"{back[q]:.5f}" for q in HARMONICS)
                        + f",{ref[i, j]:.2f},{split[i, j]:.2f}\n")


def plot_waveforms(spectra, region, title, filename):
    """B vs time per current, and its harmonic amplitudes up to the 9th."""
    fig, (ax_t, ax_f) = plt.subplots(2, 1, figsize=(8, 7))
    n_harm = max(HARMONICS) + 1
    width = FREQ / (len(spectra) + 1)
    for k, (current_amp, s) in enumerate(spectra.items()):
        times, b, amp, _ = s[region]
        ax_t.plot(times * 1e3, b, color=f"C{k}", label=f"{current_amp:g} A")
        orders = np.arange(1, n_harm)
        ax_f.bar(orders * FREQ + (k - len(spectra) / 2) * width, amp[1:n_harm],
                 width=width, color=f"C{k}")
    ax_t.set_xlabel("Time [ms]")
    ax_t.set_ylabel("B [T]")
    ax_t.set_title(f"{title} at 1000 RPM")
    ax_t.grid(True)
    ax_t.legend(fontsize="small", ncol=2)
    ax_f.set_xticks(np.arange(1, n_harm) * FREQ)
    ax_f.set_xticklabels([f"{q}\n{q * FREQ:.0f}" for q in range(1, n_harm)])
    ax_f.set_xlabel("Harmonic order / frequency [Hz]")
    ax_f.set_ylabel("B amplitude [T]")
    ax_f.grid(True, axis="y")
    fig.tight_layout()
    fig.savefig(PLOT_PREFIX + filename)


def plot_reconstruction(spectra, current_amp=0, region="tooth"):
    """B waveform vs its sum of the HARMONICS only (reference figure 2)."""
    times, b, amp, phase = spectra[current_amp][region]
    t = np.linspace(0, 1 / FREQ, 500)
    rebuilt = sum(amp[q] * np.cos(2 * np.pi * q * FREQ * t + phase[q]) for q in HARMONICS)
    plt.figure()
    plt.plot(times * 1e3, b, "b-", linewidth=2, label="FEMM")
    plt.plot(t * 1e3, rebuilt, "r-", linewidth=2,
             label="harmonics " + ", ".join(map(str, HARMONICS)))
    plt.xlabel("Time [ms]")
    plt.ylabel("B [T]")
    plt.title(f"{region.capitalize()} B at {current_amp:g} A: FEMM vs harmonics")
    plt.grid(True)
    plt.legend()
    plt.savefig(PLOT_PREFIX + f"reconstruction_{region}_{current_amp:g}A.png")


def plot_loss_map(currents, table, title, filename):
    plt.figure(figsize=(8, 6))
    filled = plt.contourf(RPMS, currents, table, levels=10, cmap="viridis")
    lines = plt.contour(RPMS, currents, table, levels=filled.levels, colors="k",
                        linewidths=0.5)
    plt.clabel(lines, fmt="%.0f", fontsize="small")
    plt.colorbar(filled, label="Iron loss [W]")
    plt.xlabel("RPM")
    plt.ylabel("Current [A]")
    plt.title(title)
    plt.savefig(PLOT_PREFIX + filename)


if __name__ == "__main__":
    spectra = harmonics(load_results())
    currents, ref, split = loss_maps(spectra)
    print_harmonics(spectra)
    print_losses(currents, ref, split)
    save_results(spectra, currents, ref, split)

    plot_waveforms(spectra, "tooth", "Stator Tooth B", "B_tooth_fft.png")
    plot_waveforms(spectra, "backiron", "Stator Back-Iron B", "B_backiron_fft.png")
    plot_reconstruction(spectra, 0, "tooth")
    plot_loss_map(currents, ref, "Iron core loss [W], tooth B over whole core (reference)",
                  "loss_map_reference.png")
    plot_loss_map(currents, split, "Iron core loss [W], teeth + back iron",
                  "loss_map_split.png")
    plt.show()
