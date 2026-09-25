import math

# ---------------------------------------
# Problem setup
# ---------------------------------------
Length = 3.300

# ---------------------------------------
# Per-region mesh sizing (model length units, inches). Each block label in
# assign_block_labels() is meshed independently at its own numeric element
# size (mi_setblockprop automesh=0) rather than all sharing one smartmesh-
# chosen size -- lets the airgap/magnet regions, which drive torque
# accuracy, be resolved much finer than the bulk steel back-iron.
# ---------------------------------------
StatorSteelMeshSize = 0.5   # was 0.05; tooth/yoke back-iron -- coarse, low field gradient
RotorSteelMeshSize = 0.5   # was 0.05; rotor lamination -- coarse, low field gradient
AirgapMeshSize = 0.5   # was 0.0025; airgap ring -- fine, resolves torque-critical field
MagnetMeshSize = 0.5   # was 0.01; magnet body -- fine
MagnetAirMeshSize = 0.5   # was 0.005; end-cap (outer) air pockets inside each pole leg -- fine
DuctMeshSize = 0.5   # was 0.05; duct/flux-barrier air pocket inside each pole leg -- coarse

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

RotorOD = 6.314
RotorID = 4.356

BridgeID = 0.370

DuctThick = 0.185
RibHeight = 0.118
RibWidth = 0.551

DistMinMag = 0.118

MagThick = 0.255
MagWidth = 0.744

Bridge = 0.056

alpha = 72.5

DuctMinDia = RotorID + 2 * BridgeID
DuctMaxDia = RotorOD - 2 * Bridge

# ---------------------------------------
# Cogging-only mesh split (dense tooth/rotor-edge bands vs. coarse
# yoke/rotor-bulk, and an unwound/Air slot body) -- opt-in via
# build_model(windings=False, split_stator_mesh=True, split_rotor_mesh=True),
# used by cogging_torque.py; every other script keeps the single-region
# defaults above.
# ---------------------------------------
SlotAirMeshSize = 0.5   # was 0.08; in, unwound slot body -- very coarse

# Radius separating the dense stator teeth from the coarse back-iron yoke.
# Must clear the tooth-bottom fillet, which bulges out to ~4.51 in with the
# geometry above (StatorID=6.375, ShoeHeight/ShoeRadius/SlotHeight as set),
# and stay under StatorOD/2=5.3 in.
ToothSplitRadius = 4.55
ToothMeshSize = 0.5   # was 0.03; in, tooth body (shoe band -> ToothSplitRadius) -- medium
YokeMeshSize = 0.5   # was 0.1; in, back-iron yoke -- very coarse

# Radius separating the thin, torque-critical tooth-shoe band facing the
# airgap from the rest of the tooth body -- just past the shoe, where the
# slot cavity (air pocket) opens up to its full SlotWidthTop. The shoe's
# outer corner (x5 in draw_stator_slot) sits at ~3.251 in; ShoeSplitMargin
# keeps the arc clear of those nodes so the mesh doesn't form slivers.
ShoeSplitMargin = 0.05       # in, beyond the shoe corner radius
ShoeSplitRadius = StatorID / 2 + ShoeHeight + ShoeRadius + ShoeSplitMargin
ShoeMeshSize = 0.5   # was 0.0025; in, tooth-shoe band -- very dense (matches airgap)

# Radius (full 0-to-SectorAngle arc) separating the dense rotor-edge band
# facing the stator from the coarse bulk rotor steel -- set to the
# mechanical bridge's own inner radius (RotorOD/2 - Bridge) so the band's
# thickness is exactly the bridge and the cut doesn't intersect the magnet
# pockets, which stay just inside DuctMaxDia/2 (~3.10 in, vs. RotorOD/2's
# 3.157 in) by design.
RotorEdgeRadius = DuctMaxDia / 2
RotorEdgeMeshSize = 0.5   # was 0.0025; in, rotor-edge band (second layer, below the skin) -- dense
RotorBulkMeshSize = 0.5   # was 0.1; in, rotor steel below RotorInnerRadius -- very coarse

# Radius of an inner rotor arc, just inside the deepest magnet-pocket point
# (~2.521 in, the magnet's inner corner) so it crosses no pocket. The layer
# between it and RotorEdgeRadius -- the steel around the magnets -- gets a
# fine mesh; only the steel below it stays at RotorBulkMeshSize.
RotorInnerRadius = 2.50
RotorMidMeshSize = 0.5   # was 0.01; in, rotor steel around the magnets -- fine

# Radius of a second, outer rotor arc splitting the edge band again: the
# thin skin between it and RotorOD (the surface directly facing the airgap)
# gets the finest rotor mesh. Must stay between RotorEdgeRadius (~3.101 in)
# and RotorOD/2 (3.157 in).
RotorSkinThick = 0.02        # in, skin thickness below RotorOD
RotorSkinRadius = RotorOD / 2 - RotorSkinThick
RotorSkinMeshSize = 0.5   # was 0.0025; in, rotor outer skin -- very dense (matches airgap)

PolePitch = 360 / Npoles        # 45 deg, rotor pole pitch
PoleHalfAngle = PolePitch / 2   # 22.5 deg

# ---------------------------------------
# Sector (periodic) reduction
# ---------------------------------------
# Smallest repeating unit: 360 / GCD(Nslots, Npoles) = 360/8 = 45 deg
# (6 stator slots + 1 rotor pole). The sector spans [0, SectorAngle],
# starting flush on the x-axis (not centered on it).
ToothPitch = 360 / Nslots                       # 7.5 deg, stator slot pitch
SectorAngle = 360 / math.gcd(Nslots, Npoles)    # 45 deg

# Sector cuts and the sliding band use anti-periodic boundaries by default,
# correct when SectorAngle spans an odd number of pole pitches (here: one
# pole, so the sector's mirror image across each cut is the opposite
# magnetic polarity). Set True only to A/B compare against a periodic
# boundary -- physically correct only if the sector spans an even number of
# pole pitches (e.g. two full poles), which is NOT the case for the default
# Nslots/Npoles here.
UsePeriodicBoundary = False

assert SectorAngle % ToothPitch == 0, \
    "sector boundary does not land on a slot-pitch multiple"
SectorSlots = int(round(SectorAngle / ToothPitch))   # 6 stator slots per sector

# The x-axis passes exactly between two slots (a tooth center), so the base
# slot unit is drawn rotated by ToothPitch/2 (3.75 deg) off the x-axis, then
# tiled by ToothPitch (7.5 deg) SectorSlots times -- 6 x 7.5 deg = 45 deg,
# landing tooth centers (not slot mouths) at the sector boundaries 0/45.
StatorBaseRotation = ToothPitch / 2   # 3.75 deg

# ---------------------------------------
# Cogging torque sweep
# ---------------------------------------
# Cogging torque repeats every 360/lcm(Nslots, Npoles) mechanical degrees --
# the angle after which the slot/pole pattern re-aligns with itself.
# (math.lcm is 3.9+; computed via gcd for compatibility with 3.8.)
NPoleSlotLCM = Nslots * Npoles // math.gcd(Nslots, Npoles)
CoggingPeriodAngle = 360 / NPoleSlotLCM   # 7.5 deg here

# Angular resolution of the sweep: how many solve steps cover one period.
CoggingStepsPerPeriod = 30
CoggingStepAngle = CoggingPeriodAngle / CoggingStepsPerPeriod

# Total mechanical angle to sweep.
CoggingSweepAngle = 7.5
CoggingNumPeriods = CoggingSweepAngle / CoggingPeriodAngle
CoggingNumSteps = int(round(CoggingNumPeriods * CoggingStepsPerPeriod))

# ---------------------------------------
# Windings / circuits (reference.txt's '19 AWG'/A-B-C setup, adapted to the
# one-pole sector: SectorSlots consecutive slots instead of all Nslots).
# ---------------------------------------
WireDiameter = 0.912        # mm -- mi_addmaterial's WireD is always mm,
                             # regardless of the document's own length units
WireConductivity = 58       # MS/m, copper (reference.txt's '19 AWG' Cduct)
WindingMeshSize = 0.5   # was 0.020; in, per-block mesh size for the coil labels

Current = 10                # A, phase current amplitude
Turns = 117                 # turns per coil side
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
                                     # reference (phdengineeringem.blogspot
                                     # .com/2018/06/toyota-prius-torque-
                                     # calculation.html).

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
# Torque-vs-current verification sweep -- reference_1809.txt's approach for
# checking this model's torque against published/measured Toyota Prius
# test data. Rotor held fixed (StepAngle=0 there, same sector-vs-full-model
# calibration as static_sweep.py's run_static_sweep()); for each current
# amplitude tested, the current's electrical PHASE (not time) is swept to
# find that amplitude's peak torque, tracing out a peak-torque-vs-current
# (Kt) curve.
# ---------------------------------------
TorqueVsCurrentAmps = [50, 75, 100, 125, 150, 200, 250]  # A, tested current levels
TorqueVsCurrentPhaseInit = 0    # deg electrical, phase sweep start
TorqueVsCurrentPhaseStep = 8    # deg electrical, phase sweep resolution
TorqueVsCurrentPhaseSteps = 22  # niterat -- sweeps Phase 0..176 deg, 23 points
