import math

MM = 1 / 25.4   # in per mm -- for dimensions taken exactly from the mm reference

# ---------------------------------------
# Problem setup_ inches
# ---------------------------------------
Length = 83.6 * MM   # 3.2913 in

# ---------------------------------------
StatorSteelMeshSize = 0.05   # tooth/yoke back-iron -- coarse, low field gradient
RotorSteelMeshSize = 0.05   # rotor lamination -- coarse, low field gradient
AirgapMeshSize = 0.0025   # airgap ring -- fine, resolves torque-critical field
MagnetMeshSize = 0.05   # magnet body -- fine
MagnetAirMeshSize = 0.15   # end-cap (outer) air pockets inside each pole leg -- fine
DuctMeshSize = 0.05   # duct/flux-barrier air pocket inside each pole leg -- coarse

# ---------------------------------------
# Stator Geometry
# ---------------------------------------
Nslots = 48

StatorOD = 10.600
StatorID = 6.375

ShoeHeight = 0.040
ShoeRadius = 0.020
SlotHeight = 1.149
SlotDia = 0.222

SlotOpen = 0.076
SlotWidthTop = 0.124
SlotWidthBot = SlotDia

StatorPitchAirgap = math.pi * StatorID / Nslots
StatorPitchSlotTop = math.pi * (StatorID + 2 * ShoeHeight) / Nslots
StatorPitchSlotBot = math.pi * (
    StatorID
    + 2 * ShoeHeight
    + 2 * ShoeRadius
    + 2 * SlotHeight
) / Nslots

# ---------------------------------------
# Rotor Geometry
# ---------------------------------------
Npoles = 8

RotorOD = 160.47 * MM   # 6.3177 in
RotorID = 111 * MM      # 4.3701 in

BridgeID = 0.370

DuctThick = 0.185
RibHeight = 0.118
RibWidth = 0.551

DistMinMag = 0.118

MagThick = 6.5 * MM    # 0.2559 in
MagWidth = 0.744

Bridge = 0.056

alpha = 72.5

DuctMinDia = RotorID + 2 * BridgeID
DuctMaxDia = RotorOD - 2 * Bridge

# Cogging-only mesh split (dense tooth/rotor-edge bands vs. coarse

SlotAirMeshSize = 0.08   # in, unwound slot body -- very coarse

ToothSplitRadius = 4.55
ToothMeshSize = 0.03   # in, tooth body (shoe band -> ToothSplitRadius) -- medium
YokeMeshSize = 0.1   # in, back-iron yoke -- very coarse


ShoeSplitMargin = 0.05       # in, beyond the shoe corner radius
ShoeSplitRadius = StatorID / 2 + ShoeHeight + ShoeRadius + ShoeSplitMargin
ShoeMeshSize = 0.0025   # in, tooth-shoe band -- very dense (matches airgap)


RotorEdgeRadius = DuctMaxDia / 2
RotorEdgeMeshSize = 0.0025   # in, rotor-edge band (second layer, below the skin) -- dense
RotorBulkMeshSize = 0.1   # in, rotor steel below RotorInnerRadius -- very coarse


RotorInnerRadius = 2.50
RotorMidMeshSize = 0.01   # in, rotor steel around the magnets -- fine


RotorSkinThick = 0.02        # in, skin thickness below RotorOD
RotorSkinRadius = RotorOD / 2 - RotorSkinThick
RotorSkinMeshSize = 0.0025   # in, rotor outer skin -- very dense (matches airgap)

PolePitch = 360 / Npoles        # 45 deg, rotor pole pitch
PoleHalfAngle = PolePitch / 2   # 22.5 deg


ToothPitch = 360 / Nslots                       # 7.5 deg, stator slot pitch
SectorAngle = 360 / math.gcd(Nslots, Npoles)    # 45 deg

UsePeriodicBoundary = False

assert SectorAngle % ToothPitch == 0, \
    "sector boundary does not land on a slot-pitch multiple"
SectorSlots = int(round(SectorAngle / ToothPitch))   # 6 stator slots per sector

StatorBaseRotation = ToothPitch / 2   # 3.75 deg

# Cogging torque sweep

NPoleSlotLCM = Nslots * Npoles // math.gcd(Nslots, Npoles)
CoggingPeriodAngle = 360 / NPoleSlotLCM   # 7.5 deg here

# Angular resolution of the sweep: how many solve steps cover one period.
CoggingStepsPerPeriod = 30
CoggingStepAngle = CoggingPeriodAngle / CoggingStepsPerPeriod

# Total mechanical angle to sweep.
CoggingSweepAngle = 7.5
CoggingNumPeriods = CoggingSweepAngle / CoggingPeriodAngle
CoggingNumSteps = int(round(CoggingNumPeriods * CoggingStepsPerPeriod))

# Material properties (steel, magnet, wire) live in materials.py.

# ---------------------------------------
# Windings / circuits (reference.txt's '19 AWG'/A-B-C setup, adapted to the
# one-pole sector: SectorSlots consecutive slots instead of all Nslots).
# ---------------------------------------
WindingMeshSize = 0.020   # in, per-block mesh size for the coil labels

Current = 10                # A, phase current amplitude
Turns = 9                 # turns per coil side, 117 for static simulations
Phase = 120                 # deg, phase A current angle at t=0 -- pairs with
                             # MaxTorqueInitialAngle below (reference_1.txt's
                             # setup) to align phase A's current with the
                             # d-axis at the start of the synchronized sweep
SpeedRPM = 1000             # RPM, for the electrical frequency below
Freq = Npoles * SpeedRPM / 120   # Hz

# One pole's winding pattern (SectorSlots consecutive slots), taken from
# reference.txt's 12-slot/2-pole 'AABBCC' repeat -- its first 6 entries.
# The sector's anti-periodic boundaries reproduce the opposite-direction
# half for the neighboring pole automatically, so only one pole needs
# listing here.
SlotCircuits = ["A", "A", "B", "B", "C", "C"]
SlotCoilDirs = [1, 1, -1, -1, 1, 1]

# ---------------------------------------
# Synchronized (max) torque sweep -- reference_1.txt's approach: rotate the
# rotor (here, via the sliding band, not mi_moverotate) at synchronous speed
# while advancing the 3-phase currents in time together, over one electrical
# period.
# ---------------------------------------
MaxTorqueNumSteps = 90               # steps/period -- matches the full
                                     # model's own 1 mech deg/step
                                     # resolution (niterat=90 in
                                     # prius_motor_full_model.py), for a
                                     # smooth torque-ripple curve
MaxTorqueStepTime = (1 / Freq) / MaxTorqueNumSteps          # s/step
MaxTorqueStepAngle = (720 / Npoles) / MaxTorqueNumSteps     # mech deg/step --
                                     # one electrical period = 720/Npoles
                                     # mechanical degrees
MaxTorqueInitialAngle = ToothPitch  # deg -- angle between two slots (360/
                                     # Nslots), aligning phase A's current
                                     # with the d-axis at the sweep's start.
                                     # Equals 7.5 deg here, matching the
                                     # initial rotor position used in the
                                     # Toyota Prius torque-calculation

# ---------------------------------------
# Max torque sweep -- static-current variant, matching the reference blog's
# "current fixed, rotor rotates" method: hold the 3-phase currents fixed at
# their MaxTorqueInitialAngle values (no time advance) and rotate ONLY the
# rotor mechanically across one electrical period, in coarse steps. Unlike
# MaxTorqueNumSteps above (which advances currents and rotor together at
# synchronous speed), this isolates the torque-angle curve at a single fixed
# current vector.
# ---------------------------------------
MaxTorqueStaticNumSteps = 15                        # steps across the sweep
MaxTorqueStaticSweepAngle = 720 / Npoles            # one electrical period,
                                                     # 90 mech deg here
MaxTorqueStaticStepAngle = MaxTorqueStaticSweepAngle / MaxTorqueStaticNumSteps

# ---------------------------------------
# Max torque sweep -- rotating-field variant: same synchronized-rotation
# calculation as MaxTorqueNumSteps above (rotor and 3-phase currents both
# advance together over one electrical period, per prius_motor_full_model.py's
# niterat/InitialAngle/StepAngle sweep), just at the coarser 15-step
# resolution the reference blog used (1 ms/step over the 15 ms period).
# ---------------------------------------
MaxTorqueRotatingFieldNumSteps = 15                 # steps/period
MaxTorqueRotatingFieldStepTime = (1 / Freq) / MaxTorqueRotatingFieldNumSteps
MaxTorqueRotatingFieldStepAngle = (720 / Npoles) / MaxTorqueRotatingFieldNumSteps

# ---------------------------------------
# TorqueVsCurrentAmps = [50, 75, 100, 125, 150, 200, 250]  # A, tested current levels
TorqueVsCurrentAmps = [250]  # A, tested current levels
TorqueVsCurrentPhaseInit = 0    # deg electrical, phase sweep start
TorqueVsCurrentPhaseStep = 8    # deg electrical, phase sweep resolution
TorqueVsCurrentPhaseSteps = 22  # niterat -- sweeps Phase 0..176 deg, 23 points
