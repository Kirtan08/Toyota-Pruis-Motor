import math

import matplotlib.pyplot as plt
import numpy as np

import materials

# Stator Geometry (mm)
Length = 83.6            # motor length
Nslots = 48
StatorOD = 269
StatorID = 161.93
ShoeHeight = 1.02        # Hs0
ShoeRadius = 0.50        # Hs1
SlotHeight = 29.16       # Hs2
SlotDia = 5.64           # 2Rs
SlotDepth = SlotDia / 2 + SlotHeight + 0.02
StatorPitchMid = math.pi * (StatorID + SlotDepth) / Nslots
SlotsOverlap = 4 + 1
HalfTurn = (Length + StatorPitchMid * SlotsOverlap + math.pi * SlotDia / 2) * 1.5

# Winding
Npoles = 8
AWG_dia = materials.WireDiameter   # mm, AWG19 diameter
Strands = 9
ParallelWires = 13
Turns = Strands * ParallelWires
PhaseLength_m = Nslots / 3 * HalfTurn / 1000 * Strands
WireResist = materials.WireResistPerMeter   # Ohm/m
PhaseResistDC_20C = WireResist * PhaseLength_m / ParallelWires

# DC resistance vs temperature
Alpha_coef = materials.CopperTempCoef   # 1/degC
Temp = np.array([20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
PhaseResistDC = PhaseResistDC_20C * (1 + Alpha_coef * (Temp - 20))

# AC resistance vs speed (skin effect)
SpeedRPM = np.array([500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000])
Freq = Npoles * SpeedRPM / 120   # Hz
Omega = 2 * math.pi * Freq
Mu0 = materials.MU0              # H/m
Sigma = materials.WireConductivity * 1e6   # S/m
Rho = 1 / Sigma                  # Ohm*m

SkinDepth = np.sqrt(2 * Rho / Omega / Mu0)   # m, one per speed
# Rows: speed, columns: temperature.
PhaseResistAC = PhaseResistDC[np.newaxis, :] * (
    1 + 1 / 9 * (SlotDepth / 1000 / SkinDepth[:, np.newaxis]) ** 2
    * (AWG_dia / 1000 / SkinDepth[:, np.newaxis]) ** 2
)


def print_results():
    print(f"HalfTurn = {HalfTurn:.4f} mm")
    print(f"Turns = {Turns}")
    print(f"PhaseLength_m = {PhaseLength_m:.4f} m")
    print(f"PhaseResistDC_20C = {PhaseResistDC_20C:.6f} Ohm")
    print(f"PhaseResistDC = {PhaseResistDC}")
    print(f"Freq = {Freq}")
    print(f"SkinDepth [mm] = {SkinDepth * 1000}")


def plot_results():
    plt.figure(1)
    for j, temp in enumerate(Temp):
        plt.plot(SpeedRPM, PhaseResistAC[:, j], ".-", label=f"{temp}")
    plt.xlabel("Speed [RPM]")
    plt.ylabel("Phase Resistance [Ohm]")
    plt.legend(title="Temp [deg C]", loc="upper left")
    plt.grid(True, which="both")
    plt.minorticks_on()
    plt.xlim(0, SpeedRPM.max())

    plt.figure(2)
    plt.plot(Freq, SkinDepth * 1000, ".-")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Skin Depth [mm]")
    plt.grid(True, which="both")
    plt.minorticks_on()
    plt.xlim(0, Freq.max())

    plt.figure(3)
    plt.plot(Temp, PhaseResistDC, ".-")
    plt.xlabel("Temperature [deg C]")
    plt.ylabel("Phase Resistance DC [Ohm]")
    plt.grid(True, which="both")
    plt.minorticks_on()
    plt.xlim(0, Temp.max())

    plt.show()


if __name__ == "__main__":
    print_results()
    plot_results()
