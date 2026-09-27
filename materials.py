"""Material definitions for the Prius motor model -- names, magnetic and
electrical properties, and the FEMM calls that add them to a document.
Import this module wherever a material property is needed so every script
uses the same values."""
import math

import femm

MU0 = 4 * math.pi * 1e-7   # H/m

# ---------------------------------------
# Material names, as used in mi_setblockprop
# ---------------------------------------
AIR = "Air"
WIRE = "19 AWG"
MAGNET = "N36Z_20"
STEEL = "M19_29G"

# ---------------------------------------
# 19 AWG copper winding wire
# ---------------------------------------
WireDiameter = 0.912        # mm -- mi_addmaterial's WireD is always mm,
                             # regardless of the document's own length units
WireConductivity = 58       # MS/m, copper at 20 deg C
WireResistPerMeter = 0.02648   # Ohm/m at 20 deg C
CopperTempCoef = 0.004041   # 1/deg C, resistance temperature coefficient
                             # (0.00393 or 0.004041)

# ---------------------------------------
# N36Z_20 magnet with temperature derating. mur stays constant while Hc (and
# so Br = mu0*mur*Hc) falls linearly with temperature. FEMM takes Hc as a
# positive magnitude; the direction comes from each block's magdir.
# ---------------------------------------
MagnetMur = 1.03
MagnetConductivity = 0.667  # MS/m
MagnetHcRef = 920000        # A/m, coercivity at MagnetRefTemp
MagnetRefTemp = 20          # deg C
MagnetTemp = 50             # deg C, operating magnet temperature
MagnetHcTempCoef = -0.5     # %/deg C, dHc/dT
MagnetHc = MagnetHcRef * (1 + MagnetHcTempCoef / 100 * (MagnetTemp - MagnetRefTemp))
# 920000 * (1 - 0.005*30) = 782000 A/m (9.83 kOe) at 50 deg C
MagnetBr = MU0 * MagnetMur * MagnetHc   # T

# ---------------------------------------
# M19_29G lamination steel
# ---------------------------------------
SteelConductivity = 1.9     # MS/m
SteelLamThickness = 0.34    # mm
# Lamination stacking factor, applied to the solid-material B-H curve below
# (flux parallel to the laminations); FEMM's own LamFill is then set to 1 so
# it isn't applied twice.
StackingFactor = 0.94

# Solid-material B-H curve (T, A/m).
SteelBPoints = [
    0, 0.05, 0.1, 0.15, 0.36, 0.54, 0.65, 0.99, 1.2, 1.28, 1.33, 1.36,
    1.44, 1.52, 1.58, 1.63, 1.67, 1.8, 1.9, 2, 2.1, 2.3, 2.5,
    2.563994494, 3.779889874,
]
SteelHPoints = [
    0, 22.28, 25.46, 31.83, 47.74, 63.66, 79.57, 159.15, 318.3, 477.46,
    636.61, 795.77, 1591.5, 3183, 4774.6, 6366.1, 7957.7, 15915, 31830,
    111407, 190984, 350135, 509252, 560177.2, 1527756,
]


def laminated_b(b_solid, h_solid, sf=None):
    """Effective flux density of a lamination stack with flux parallel to
    the laminations: Bparallel = SF*Bsolid + (1 - SF)*mu0*Hsolid. H is the
    same in the iron and the gaps between laminations, so it is unchanged."""
    if sf is None:
        sf = StackingFactor
    return sf * b_solid + (1 - sf) * MU0 * h_solid


def setup_materials():
    """Add every material above to the open FEMM document."""
    femm.mi_getmaterial(AIR)

    femm.mi_addmaterial(
        WIRE, 1.0, 1.0, 0, 0, WireConductivity, 0, 0, 0, 3, 0, 0, 1,
        WireDiameter,
    )

    femm.mi_addmaterial(MAGNET, MagnetMur, MagnetMur, MagnetHc, 0, MagnetConductivity)

    # LamFill = 1: the stacking factor is applied to the B-H points instead,
    # so FEMM must not apply it a second time.
    femm.mi_addmaterial(STEEL, 0, 0, 0, 0, SteelConductivity, SteelLamThickness, 0, 1, 0)
    for b, h in zip(SteelBPoints, SteelHPoints):
        femm.mi_addbhpoint(STEEL, laminated_b(b, h), h)


if __name__ == "__main__":
    print(f"Magnet at {MagnetTemp} deg C: Hc = {MagnetHc:.0f} A/m "
          f"({MagnetHc * 4 * math.pi / 1e6:.2f} kOe), Br = {MagnetBr:.3f} T")
