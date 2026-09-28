# Toyota Prius 2004 Motor: FEMM Sector Model

A 2D magnetostatic finite-element model of the 2004 Toyota Prius traction motor
(interior permanent magnet, V-shaped magnets, 48 slots / 8 poles), built and
solved in [FEMM](https://www.femm.info/) through its Python API (`pyfemm`).

To save solve time, the model simulates one **45° sector** (6 slots, 1 pole)
instead of the full machine. Anti-periodic boundaries on the sector cuts stand
in for the rest of the motor, and a **sliding band** in the airgap rotates the
rotor without redrawing the geometry.

## Requirements

- Windows with [FEMM 4.2](https://www.femm.info/wiki/Download) installed
- Python 3.8+
- Python packages: `pyfemm`, `matplotlib`

```
pip install pyfemm matplotlib
```

## Repository layout

| File | Purpose |
| --- | --- |
| [`config.py`](config.py) | All parameters: geometry, mesh sizes, windings, and sweep settings |
| [`materials.py`](materials.py) | Material names and properties (M19_29G B-H curve, N36Z_20 magnet with temperature derating, 19 AWG copper) and `setup_materials()` |
| [`simulation.py`](simulation.py) | Builds the model: geometry, windings, boundary conditions. Also runs on its own to build, save and solve one model (`ToyotaPrius.FEM`) |
| [`cogging_torque.py`](cogging_torque.py) | Cogging torque sweep (no current) |
| [`static_current_rotor_spinning.py`](static_current_rotor_spinning.py) | Torque vs rotor position with fixed phase currents |
| [`Locked_rotor_test.py`](Locked_rotor_test.py) | Peak torque vs current amplitude (Kt curve) |
| [`ironLoss_P2.py`](ironLoss_P2.py) | Rotating-rotor sweep at MTPA: stator tooth / back-iron B waveforms, torque, flux linkage and phase voltages |
| [`ironLoss_P3.py`](ironLoss_P3.py) | FFT of the P2 B waveforms and iron loss map vs rpm and current (no FEMM) |
| [`flux_map.py`](flux_map.py) | ψd/ψq lookup tables over an Id × Iq grid, and apparent vs incremental inductance |
| [`efficiency_map.py`](efficiency_map.py) | Copper + iron loss, output/input power and efficiency maps vs speed and torque (no FEMM) |
| [`notes/`](notes/) | Dated working notes |

## Running a simulation

Each script opens FEMM, builds the model, then **pauses so you can inspect the
geometry**. Press Enter in the terminal to start the sweep. Results are saved as
CSV and PNG files, and the plots are shown when the run finishes.

### Cogging torque

```
python cogging_torque.py
```

Runs with no current in the windings (slots filled with air) and rotates the
rotor over one cogging period (7.5° mechanical, 30 steps). The stator and rotor
steel are split into extra regions so the areas facing the airgap can be meshed
more finely than the bulk steel.

Output goes to `Cogging_outputs/`: the model (`.FEM`), results (`.csv`), plot
(`.png`), and an image of the initial mesh (`.bmp`).

![Cogging torque](Cogging_outputs/Cogging_outputs.png)

### Torque vs rotor position (fixed currents)

```
python static_current_rotor_spinning.py
```

Holds the phase currents fixed (`I_A`, `I_B`, `I_C` at the top of the script)
and rotates only the rotor across one electrical period (90° mechanical, 15
steps). Prints the maximum torque and the angle where it occurs.

Output: `static_current_rotor_spinning.csv` and
`static_current_rotor_spinning_vs_angle.png`.

![Torque vs rotor position](static_current_rotor_spinning_vs_angle.png)

### Torque vs current (Kt curve)

```
python Locked_rotor_test.py          # all currents in config.TorqueVsCurrentAmps
python Locked_rotor_test.py 150      # a single current amplitude, in A
```

Holds the rotor fixed and, for each current amplitude, sweeps the current phase
from 0° to 176° electrical to find the peak torque. This gives the peak-torque
vs current curve, for comparison with published Prius data.

Output: `torque_vs_current.csv` (every point), `torque_vs_current_kt.csv`
(peaks), and the plots `torque_vs_current_phase.png` and
`torque_vs_current_kt.png`.

### Iron loss sweep (rotating rotor, MTPA)

```
python ironLoss_P2.py          # run the FEMM sweep, then plot
python ironLoss_P2.py --plot   # replot the saved CSV only
```

Turns the rotor at 1000 rpm, synchronized with the 3-phase currents, over one
electrical cycle (90 steps) for each current from 0 to 250 A at its MTPA
current angle. Each step records the phase currents, flux linkages and
voltages, the torque, B in the middle of a stator tooth and of the back iron
(the waveforms for the iron loss), and the magnet operating point.

Phase voltages are plotted three ways:

- `ironLoss_P2_emf.png`: induced voltage (back-EMF) `dψ/dt`, from the flux
  linkage.
- `ironLoss_P2_voltage_resistive.png`: FEMM's circuit voltage. In a
  magnetostatic solve this is only the resistive drop R·i over the stack
  length, not an induced voltage. The reference Octave script plots this
  quantity as "Induced Phase Voltages".
- `ironLoss_P2_voltage_terminal.png`: terminal voltage R·i + dψ/dt.

FEMM's R is about 13× too high because the model treats each turn as one
wire, while the real coil has 13 wires in parallel. So the R·i values are
overstated. The back-EMF is not affected. See the
[2026-09-28 notes](notes/Toyota_prius_2004_2026-09-26.md#update-2026-09-28).

Output: `ironLoss_P2.csv` (every step), `ironLoss_P2_kt.csv` (per-current
summary), and the `ironLoss_P2_*.png` plots.

### Iron loss map (FFT of the stator B)

```
python ironLoss_P3.py
```

No FEMM run: reads `ironLoss_P2.csv`, splits the tooth and back-iron B
waveforms into harmonics with an FFT, and applies the loss model fitted in
`ironLoss_P1.py` to the 1st, 3rd, 5th, 7th and 9th harmonics. The waveforms are
assumed not to change with speed, so the loss is mapped from 200 to 6000 rpm
by scaling the frequency. Two maps:

- `ironLoss_P3_loss_map_reference.png`: tooth B applied to the whole stator
  core, as in the reference study.
- `ironLoss_P3_loss_map_split.png`: tooth B over the teeth and back-iron B over
  the back iron.

Output: `ironLoss_P3.csv` (harmonic amplitudes and both losses per current
and rpm), the FFT plots `ironLoss_P3_B_tooth_fft.png` and
`ironLoss_P3_B_backiron_fft.png`, and `ironLoss_P3_reconstruction_tooth_0A.png`
(the 0 A waveform rebuilt from its harmonics).

### Efficiency map

```
python efficiency_map.py
```

No FEMM run: combines the mean torque from `ironLoss_P2_kt.csv`, the iron
loss from `ironLoss_P3.csv` (teeth + back iron by default, `IRON_LOSS_COLUMN`),
and the copper loss 3/2·R_ac·I_peak² with R_ac from
`phase_resistance_calculation.phase_resistance_ac()` at 140 °C. Efficiency is
P_out / (P_out + P_copper + P_iron). Points with P_in above 60 kW are blanked.

Each map shows the 12.5 / 25 / 37.5 / 50 kW output power lines and a red
**voltage limit at MTPA** line: the speed where ω·ψ_peak + R·I reaches
500 V / √3 (the Prius's maximum DC link voltage). The model has no field
weakening, so points to the right of that line can't be reached at MTPA and
their losses are not valid.

Output: `efficiency_map.csv` and `efficiency_map_{copper_loss, iron_loss,
total_loss, output_power, input_power, efficiency}.png`.

### Flux linkage lookup tables (Id, Iq)

```
python flux_map.py          # FEMM sweep, 182 solves (~32 min), then tables and plots
python flux_map.py --plot   # rebuild from flux_map.csv
```

Solves a grid of Id (−600 to +50 A) and Iq (0 to 600 A) in 50 A steps, with
the magnets on and the rotor at the locked-rotor position, and transforms the
phase flux linkages to ψd(Id, Iq), ψq(Id, Iq). From these it derives:

- **apparent** inductances: Ld = (ψd − ψm)/Id, Lq = ψq/Iq. They give the flux
  itself, so they are the ones for torque and steady-state voltage.
- **incremental** inductances: ∂ψd/∂Id, ∂ψq/∂Iq and the cross terms
  ∂ψd/∂Iq, ∂ψq/∂Id. They give how the flux changes, so they are the ones for
  di/dt, current ripple and current-controller gains.

Output:
- `flux_map_lut_<name>.csv`: lookup tables as matrices (Id down, Iq across)
  for ψd, ψq, torque and each inductance;
- `flux_map_surface_<name>.png`: 3D surfaces of the same tables;
- `flux_map_{tables, psi_curves, apparent_vs_incremental, cross_coupling,
  mtpa_comparison}.png`.

`ROTOR_POSITIONS > 1` averages each point over one slot pitch to remove the
slotting ripple. A finer 25 A grid below 100 A is commented out in the script.

## Model details

- **Units:** inches; stack length 3.3 in.
- **Materials:** M19_29G electrical steel (nonlinear B-H curve), N36Z_20
  magnets, 19 AWG copper windings.
- **Windings:** 3-phase, 117 turns per coil side, one pole of the `AABBCC`
  slot pattern.
- **Torque:** calculated with FEMM's air-gap integral on the sliding band
  (`mo_gapintegral`).
- **Mesh:** every region has its own mesh size in `config.py`. Smartmesh
  (automatic sizing) was found to badly underestimate torque, so the airgap
  mesh is sized manually.

> **Note:** all mesh sizes in `config.py` are currently set to a coarse `0.5` in.
> The comment next to each one gives its previous, finer value (for example the
> airgap was `0.0025`). Restore those values for accurate results. Coarse meshes
> run quickly but are only suitable for checking that the model works.

## Status

The cogging torque result is about 35% off published results from commercial
software, even after mesh refinement. See
[`notes/Toyota_prius_2004_2026-09-26.md`](notes/Toyota_prius_2004_2026-09-26.md).


