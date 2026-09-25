import femm
import math
import config

# ---------------------------------------
# Helper functions for degree trig
# ---------------------------------------
def sind(x):
    return math.sin(math.radians(x))

def cosd(x):
    return math.cos(math.radians(x))

def tand(x):
    return math.tan(math.radians(x))


# ---------------------------------------
# Geometry helper functions
# ---------------------------------------
def arc(radius, angle_start, angle_end, cx=0, cy=0, group=None):
    """Draw an arc of the given radius from angle_start to angle_end (deg)."""
    x_start = cx + radius * cosd(angle_start)
    y_start = cy + radius * sind(angle_start)
    x_end = cx + radius * cosd(angle_end)
    y_end = cy + radius * sind(angle_end)

    femm.mi_addnode(x_start, y_start)
    femm.mi_addnode(x_end, y_end)
    femm.mi_addarc(x_start, y_start, x_end, y_end, angle_end - angle_start, 1)

    if group is not None:
        femm.mi_selectarcsegment((x_start + x_end) / 2, (y_start + y_end) / 2)
        femm.mi_setgroup(group)
        femm.mi_clearselected()

    return (x_start, y_start), (x_end, y_end)


def radial_line(angle, r_start, r_end, cx=0, cy=0, group=None):
    """Draw a straight radial cut at the given angle (deg) from r_start to r_end."""
    x_start = cx + r_start * cosd(angle)
    y_start = cy + r_start * sind(angle)
    x_end = cx + r_end * cosd(angle)
    y_end = cy + r_end * sind(angle)

    femm.mi_addnode(x_start, y_start)
    femm.mi_addnode(x_end, y_end)
    femm.mi_addsegment(x_start, y_start, x_end, y_end)

    if group is not None:
        femm.mi_selectsegment((x_start + x_end) / 2, (y_start + y_end) / 2)
        femm.mi_setgroup(group)
        femm.mi_clearselected()

    return (x_start, y_start), (x_end, y_end)


def mirror_about_line(x, y, ax, ay):
    """Reflect (x, y) about the line through the origin and (ax, ay) --
    vector reflection, taking the same axis POINT mi_mirror2 is given
    (rather than re-deriving its angle independently), so a manually
    computed mirror image lands on the exact same node FEMM's own mirror
    creates, instead of a numerically-adjacent-but-distinct one (which
    meshes into a degenerate sliver element and can crash the solver)."""
    d2 = ax * ax + ay * ay
    dot = x * ax + y * ay
    return (2 * dot * ax / d2 - x, 2 * dot * ay / d2 - y)


def rotate_point(x, y, angle):
    """Rotate a point (x, y) by angle (deg) about the origin."""
    return (
        x * cosd(angle) - y * sind(angle),
        x * sind(angle) + y * cosd(angle),
    )


# ---------------------------------------
# FEMM group numbers. Every segment/arc making up a stator slot is tagged
# with SLOT_GROUP via mi_setgroup, so mi_selectgroup(SLOT_GROUP) selects all
# of them at once regardless of which slot (or how many) drew them.
# ---------------------------------------
SLOT_GROUP = 1
MAGNET_GROUP = 2
ROTOR_STEEL_GROUP = 3
STATOR_STEEL_GROUP = 4
AIRGAP_GROUP = 5
BAND_GROUP = 6

# Name of the anti-periodic sliding-band boundary used to rotate the rotor
# at solve time (via mi_modifyboundprop) without redrawing geometry.
SLIDING_BAND_NAME = "mySlidingBand"

# Sliding-band interface radius: midway across the physical airgap between
# the rotor OD and the stator ID.
SLIDING_BAND_RADIUS = (config.RotorOD / 2 + config.StatorID / 2) / 2


# ---------------------------------------
# Sector geometry (outline only -- no magnets, BCs, or materials yet)
# ---------------------------------------
def draw_sector_outline():
    # Sector spans [0, SectorAngle], starting flush on the x-axis.
    # Stator ID is NOT drawn as a smooth arc: the slot mouths tiled by
    # draw_stator_slots() form that boundary instead.
    arc(config.StatorOD / 2, 0, config.SectorAngle)
    arc(config.RotorOD / 2, 0, config.SectorAngle)
    arc(config.RotorID / 2, 0, config.SectorAngle)
    # Sliding-band arc at mid-airgap radius -- lets the rotor be rotated at
    # solve time (via mi_modifyboundprop) without redrawing geometry.
    arc(SLIDING_BAND_RADIUS, 0, config.SectorAngle, group=BAND_GROUP)

    for angle in (0, config.SectorAngle):
        radial_line(angle, config.StatorID / 2, config.StatorOD / 2)
        radial_line(angle, config.RotorID / 2, config.RotorOD / 2)
        # Connect the sliding-band arc to the stator ID so the band's
        # endpoints aren't left dangling. The rotor side (RotorOD -> band)
        # is deliberately left angularly open and unmeshed: per FEMM's
        # sliding-band technique (femm.info/wiki/SlidingBand), that thin
        # strip is NOT a normal triangulated block region -- FEMM builds its
        # own internal quadrilateral mesh there by interpolating between the
        # rotor-side and stator-side node rings, using the (anti)periodic
        # air-gap boundary condition applied to the arcs on both sides of
        # the strip (here: RotorOD and mySlidingBand). Giving it a block
        # label/radial closure, as a previous attempt here did, fights that
        # mechanism rather than fixing anything.
        radial_line(angle, SLIDING_BAND_RADIUS, config.StatorID / 2, group=AIRGAP_GROUP)


def tooth_centerline_point(angle_offset=0):
    """The point where a slot's half tooth-bottom fillet ends -- i.e. the
    tooth centerline at ToothPitch/2 off the slot center -- used as the
    second point (with the origin) defining the mirror line between
    neighboring half-teeth. Independent of Nslots/SectorSlots."""
    c = config

    x7 = c.StatorID / 2 + c.ShoeHeight + c.ShoeRadius + c.SlotHeight
    y7 = (c.StatorPitchSlotBot - c.SlotWidthBot) / 2
    x8, y8 = x7, -y7
    x9, y9 = rotate_point(x8, y8, c.ToothPitch)
    cx, cy = (x7 + x9) / 2, (y7 + y9) / 2
    vx, vy = x7 - cx, y7 - cy
    xCutR = cx + (vx * cosd(90) - vy * sind(90))
    yCutR = cy + (vx * sind(90) + vy * cosd(90))

    return rotate_point(xCutR, yCutR, angle_offset)


# ---------------------------------------
# Stator slot (single, self-contained slot unit -- computed centered on
# theta=0 with the x-axis splitting it symmetrically top/bottom, then
# rotated in-place by angle_offset before being drawn -- no FEMM copy/rotate
# needed, and no reaching into any neighboring slot's geometry).
# ---------------------------------------
def draw_stator_slot(angle_offset=0):
    c = config

    x1 = (c.StatorID / 2) * math.cos(
        math.asin((c.StatorPitchAirgap - c.SlotOpen) / c.StatorID)
    )
    y1 = (c.StatorPitchAirgap - c.SlotOpen) / 2
    x2, y2 = x1, -y1

    x3 = c.StatorID / 2 + c.ShoeHeight
    y3 = (c.StatorPitchAirgap - c.SlotOpen) / 2

    x4 = c.StatorID / 2 + c.ShoeHeight
    y4 = -(c.StatorPitchAirgap - c.SlotOpen) / 2

    x5 = c.StatorID / 2 + c.ShoeHeight + c.ShoeRadius
    y5 = (c.StatorPitchSlotTop - c.SlotWidthTop) / 2

    x6 = c.StatorID / 2 + c.ShoeHeight + c.ShoeRadius
    y6 = -(c.StatorPitchSlotTop - c.SlotWidthTop) / 2

    x7 = c.StatorID / 2 + c.ShoeHeight + c.ShoeRadius + c.SlotHeight
    y7 = (c.StatorPitchSlotBot - c.SlotWidthBot) / 2

    x8, y8 = x7, -y7

    # x7 does NOT close directly to x8: the deep corner is rounded off by the
    # TOOTH's own bottom fillet, which spans a full ToothPitch across to the
    # matching corner of the NEXT slot (rotate x8 by +ToothPitch to get that
    # far corner, x9) -- confirmed against the working reference model
    # (chevy_bolt_v1.py), where the same p7->p9 (via Gamma rotation) pattern
    # is used. A standalone slot with no real neighbor to reach has no such
    # far corner, so it closes only to the HALF of that fillet nearest to it:
    # the arc's own true angular midpoint, which lands exactly on the tooth's
    # centerline (ToothPitch/2 off this slot's center).
    xCutR, yCutR = tooth_centerline_point()
    # Mirror (negate y) for the identical half-fillet closure on the left.
    xCutL, yCutL = xCutR, -yCutR

    # Model the slot with the x-axis (angle 0) as its own center, then rotate
    # every point into place by angle_offset before drawing anything.
    (x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5) = (
        rotate_point(x1, y1, angle_offset),
        rotate_point(x2, y2, angle_offset),
        rotate_point(x3, y3, angle_offset),
        rotate_point(x4, y4, angle_offset),
        rotate_point(x5, y5, angle_offset),
    )
    (x6, y6), (x7, y7), (x8, y8) = (
        rotate_point(x6, y6, angle_offset),
        rotate_point(x7, y7, angle_offset),
        rotate_point(x8, y8, angle_offset),
    )
    (xCutR, yCutR), (xCutL, yCutL) = (
        rotate_point(xCutR, yCutR, angle_offset),
        rotate_point(xCutL, yCutL, angle_offset),
    )

    for x, y in ((x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5),
                 (x6, y6), (x7, y7), (x8, y8), (xCutR, yCutR), (xCutL, yCutL)):
        femm.mi_addnode(x, y)

    arc_angle = (
        2 * math.asin((c.StatorPitchAirgap - c.SlotOpen) / c.StatorID)
        * 180 / math.pi
    )
    femm.mi_addarc(x2, y2, x1, y1, arc_angle, 1)
    femm.mi_addsegment(x1, y1, x3, y3)
    femm.mi_addsegment(x2, y2, x4, y4)
    femm.mi_addarc(x3, y3, x5, y5, 90, 1)
    femm.mi_addarc(x6, y6, x4, y4, 90, 1)
    femm.mi_addsegment(x5, y5, x7, y7)
    femm.mi_addsegment(x6, y6, x8, y8)
    # Half of the tooth-bottom fillet on each side, ending at the tooth's own
    # centerline -- self-contained, no reach into a neighboring slot.
    femm.mi_addarc(x7, y7, xCutR, yCutR, 90, 1)
    femm.mi_addarc(xCutL, yCutL, x8, y8, 90, 1)

    femm.mi_selectarcsegment((x1 + x2) / 2, (y1 + y2) / 2)
    femm.mi_selectsegment((x1 + x3) / 2, (y1 + y3) / 2)
    femm.mi_selectsegment((x2 + x4) / 2, (y2 + y4) / 2)
    femm.mi_selectarcsegment((x3 + x5) / 2, (y3 + y5) / 2)
    femm.mi_selectarcsegment((x4 + x6) / 2, (y4 + y6) / 2)
    femm.mi_selectsegment((x5 + x7) / 2, (y5 + y7) / 2)
    femm.mi_selectsegment((x6 + x8) / 2, (y6 + y8) / 2)
    femm.mi_selectarcsegment((x7 + xCutR) / 2, (y7 + yCutR) / 2)
    femm.mi_selectarcsegment((xCutL + x8) / 2, (yCutL + y8) / 2)
    femm.mi_setgroup(SLOT_GROUP)
    femm.mi_clearselected()


def draw_stator_half_unit(angle_offset=0):
    """Debug helper: same points/edges as draw_stator_slot, but only the
    upper (y >= 0) half chain -- x1, x3, x5, x7, and the half tooth-fillet
    closure at xCutR. Nothing below the x-axis is drawn."""
    c = config

    x0 = c.StatorID / 2
    y0 = 0

    x1 = (c.StatorID / 2) * math.cos(
        math.asin((c.StatorPitchAirgap - c.SlotOpen) / c.StatorID)
    )
    y1 = (c.StatorPitchAirgap - c.SlotOpen) / 2

    x3 = c.StatorID / 2 + c.ShoeHeight
    y3 = (c.StatorPitchAirgap - c.SlotOpen) / 2

    x5 = c.StatorID / 2 + c.ShoeHeight + c.ShoeRadius
    y5 = (c.StatorPitchSlotTop - c.SlotWidthTop) / 2

    x7 = c.StatorID / 2 + c.ShoeHeight + c.ShoeRadius + c.SlotHeight
    y7 = (c.StatorPitchSlotBot - c.SlotWidthBot) / 2

    xCutR, yCutR = tooth_centerline_point()

    arc_angle_half = (
        math.asin((c.StatorPitchAirgap - c.SlotOpen) / c.StatorID)
        * 180 / math.pi
    )

    (x0, y0), (x1, y1), (x3, y3), (x5, y5), (x7, y7), (xCutR, yCutR) = (
        rotate_point(x0, y0, angle_offset),
        rotate_point(x1, y1, angle_offset),
        rotate_point(x3, y3, angle_offset),
        rotate_point(x5, y5, angle_offset),
        rotate_point(x7, y7, angle_offset),
        rotate_point(xCutR, yCutR, angle_offset),
    )

    for x, y in ((x0, y0), (x1, y1), (x3, y3), (x5, y5), (x7, y7), (xCutR, yCutR)):
        femm.mi_addnode(x, y)

    # Half of the airgap mouth arc, from the x-axis (tooth centerline on the
    # stator ID) to x1 -- the counterpart of the segments/arcs already drawn
    # for the rest of the chain, which was missing before.
    femm.mi_addarc(x0, y0, x1, y1, arc_angle_half, 1)
    femm.mi_addsegment(x1, y1, x3, y3)
    femm.mi_addarc(x3, y3, x5, y5, 90, 1)
    femm.mi_addsegment(x5, y5, x7, y7)
    femm.mi_addarc(x7, y7, xCutR, yCutR, 90, 1)

    femm.mi_selectarcsegment((x0 + x1) / 2, (y0 + y1) / 2)
    femm.mi_selectsegment((x1 + x3) / 2, (y1 + y3) / 2)
    femm.mi_selectarcsegment((x3 + x5) / 2, (y3 + y5) / 2)
    femm.mi_selectsegment((x5 + x7) / 2, (y5 + y7) / 2)
    femm.mi_selectarcsegment((x7 + xCutR) / 2, (y7 + yCutR) / 2)
    femm.mi_setgroup(SLOT_GROUP)
    femm.mi_clearselected()


# ---------------------------------------
# Rotor pole (single V-magnet leg -- computed centered on theta=0, upper
# side only, matching the reference script's x11-x18 chain, then rotated
# in-place by angle_offset before being drawn).
# ---------------------------------------
def _rotor_pole_leg_points(angle_offset=0, mirror=False):
    """Compute x11..x18 (see draw_rotor_pole_leg) centered on theta=0,
    mirrored (y negated) if requested, then rotated into place. Returns
    {11: (x, y), ..., 18: (x, y)}. Factored out of draw_rotor_pole_leg so
    label-placement code can reuse the exact same points."""
    c = config

    x11 = c.DuctMinDia / 2
    y11 = 0

    x12 = x11 + c.DuctThick / sind(c.alpha)
    y12 = 0

    x13 = x12 + c.DistMinMag / 2 / tand(c.alpha)
    y13 = c.DistMinMag / 2

    x14 = x13 - c.DuctThick * cosd(90 - c.alpha)
    y14 = y13 + c.DuctThick * sind(90 - c.alpha)

    x15 = x14 + c.MagWidth * cosd(c.alpha)
    y15 = y14 + c.MagWidth * sind(c.alpha)

    x16 = x13 + c.MagWidth * cosd(c.alpha)
    y16 = y13 + c.MagWidth * sind(c.alpha)

    x17 = x13 - c.MagThick * cosd(90 - c.alpha)
    y17 = y13 + c.MagThick * sind(90 - c.alpha)

    x18 = x16 - c.MagThick * cosd(90 - c.alpha)
    y18 = y16 + c.MagThick * sind(90 - c.alpha)

    pts = {
        11: (x11, y11), 12: (x12, y12), 13: (x13, y13), 14: (x14, y14),
        15: (x15, y15), 16: (x16, y16), 17: (x17, y17), 18: (x18, y18),
    }
    if mirror:
        pts = {k: (x, -y) for k, (x, y) in pts.items()}
    return {k: rotate_point(x, y, angle_offset) for k, (x, y) in pts.items()}


def rotor_pole_leg_label_points(angle_offset=0, mirror=False):
    """Interior points for the 3 closed regions inside one V-magnet leg:
    the duct/flux-barrier air pocket (quad 11-12-13-14), the magnet
    rectangle (13-16-18-17; 14 and 15 lie on its edges, not corners), and
    the semicircular end-cap air pocket (bigon between chord 15-16 and its
    180 deg arc)."""
    p = _rotor_pole_leg_points(angle_offset, mirror)

    duct_pt = (
        sum(p[k][0] for k in (11, 12, 13, 14)) / 4,
        sum(p[k][1] for k in (11, 12, 13, 14)) / 4,
    )
    magnet_pt = (
        sum(p[k][0] for k in (13, 16, 18, 17)) / 4,
        sum(p[k][1] for k in (13, 16, 18, 17)) / 4,
    )

    # mi_addarc's start/end order differs by `mirror` (matches the actual
    # calls in draw_rotor_pole_leg below); the peak of a 180 deg CCW arc
    # from S to E is at center + (dy/2, -dx/2), so a point halfway between
    # the chord midpoint and the peak is safely interior to the half-disk.
    if not mirror:
        (sx, sy), (ex, ey) = p[16], p[15]
    else:
        (sx, sy), (ex, ey) = p[15], p[16]
    cx, cy = (sx + ex) / 2, (sy + ey) / 2
    dx, dy = ex - sx, ey - sy
    endcap_pt = (cx + dy / 4, cy - dx / 4)

    return duct_pt, magnet_pt, endcap_pt


def rotor_pole_leg_magdirection(angle_offset=0, mirror=False):
    """Magnetization angle (deg, absolute, 0 = +x axis) for the magnet
    region: perpendicular to the magnet's long axis, pointing radially
    outward for mirror=False and radially inward for mirror=True, forming
    a coherent V-pole with opposite poles facing each other."""
    # Invert the magnetization direction between the two pole legs so the
    # pair is flipped relative to the original orientation.
    if mirror:
        local_angle = config.alpha - 55
    else:
        local_angle = config.alpha + 90
    return (local_angle + angle_offset) % 360


def draw_rotor_pole_leg(angle_offset=0, mirror=False):
    p = _rotor_pole_leg_points(angle_offset, mirror)
    (x11, y11), (x12, y12), (x13, y13), (x14, y14) = p[11], p[12], p[13], p[14]
    (x15, y15), (x16, y16), (x17, y17), (x18, y18) = p[15], p[16], p[17], p[18]

    for x, y in ((x11, y11), (x12, y12), (x13, y13), (x14, y14),
                 (x15, y15), (x16, y16), (x17, y17), (x18, y18)):
        femm.mi_addnode(x, y)

    femm.mi_addsegment(x12, y12, x13, y13)
    femm.mi_addsegment(x11, y11, x14, y14)
    femm.mi_addsegment(x13, y13, x14, y14)
    femm.mi_addsegment(x13, y13, x16, y16)
    femm.mi_addsegment(x15, y15, x16, y16)
    if mirror:
        femm.mi_addarc(x15, y15, x16, y16, 180, 1)
    else:
        femm.mi_addarc(x16, y16, x15, y15, 180, 1)
    femm.mi_addsegment(x14, y14, x17, y17)
    femm.mi_addsegment(x18, y18, x15, y15)
    femm.mi_addsegment(x18, y18, x17, y17)

    femm.mi_selectsegment((x11 + x14) / 2, (y11 + y14) / 2)
    femm.mi_selectsegment((x12 + x13) / 2, (y12 + y13) / 2)
    femm.mi_selectsegment((x13 + x14) / 2, (y13 + y14) / 2)
    femm.mi_selectsegment((x14 + x17) / 2, (y14 + y17) / 2)
    femm.mi_selectsegment((x13 + x16) / 2, (y13 + y16) / 2)
    femm.mi_selectsegment((x15 + x16) / 2, (y15 + y16) / 2)
    femm.mi_selectsegment((x15 + x18) / 2, (y15 + y18) / 2)
    femm.mi_selectarcsegment((x15 + x16) / 2, (y15 + y16) / 2)
    femm.mi_selectsegment((x17 + x18) / 2, (y17 + y18) / 2)
    femm.mi_setgroup(MAGNET_GROUP)
    femm.mi_clearselected()


def draw_stator_slots():
    # 6 slot centers at StatorBaseRotation + k*ToothPitch (k=0..5), landing
    # exactly on tooth centers at the sector boundaries (0 and 45 deg).
    for k in range(config.SectorSlots):
        draw_stator_slot(config.StatorBaseRotation + k * config.ToothPitch)


def select_all_slot_segments():
    """Select every segment/arc tagged SLOT_GROUP, in one call, regardless of
    how many slots drew them."""
    femm.mi_selectgroup(SLOT_GROUP)


# ---------------------------------------
# Cogging-only mesh region splits: mesh-only cuts (no boundary/material
# change) that break the single stator-steel/rotor-steel region into a
# dense sub-region facing the airgap and a coarse bulk sub-region, so each
# can be meshed independently. Both halves stay M19_29G either way -- these
# functions must run AFTER the slots/magnets are drawn (so the split radius/
# angles can be checked against the fillet and magnet-pocket bulges) and
# BEFORE setup_boundary_conditions (which must re-pair whatever anti-
# periodic segments/arcs these cuts split).
# ---------------------------------------
def draw_stator_tooth_yoke_split(radius):
    """Full-sector arc at `radius`, separating the dense tooth region
    (below) from the coarse back-iron yoke (above). `radius` must clear the
    tooth-bottom fillet bulge and stay under StatorOD/2."""
    arc(radius, 0, config.SectorAngle, group=STATOR_STEEL_GROUP)


def draw_stator_shoe_split(radius):
    """Full-sector arc at `radius`, just past the tooth shoes, separating
    the dense shoe band facing the airgap (below) from the medium tooth body
    (above). It crosses every slot cavity, cutting a thin air strip off
    the bottom of each slot body; below the arc each tooth tip is its own
    isolated region (the slot mouths separate them). Both are labelled per
    tooth/slot in assign_block_labels."""
    arc(radius, 0, config.SectorAngle, group=STATOR_STEEL_GROUP)


def draw_rotor_edge_split(radius):
    """Full-sector arc at `radius`, separating the dense rotor-edge band
    facing the stator (above) from the coarse bulk rotor steel (below).
    `radius` must clear the magnet-pocket bulge (both legs) and stay under
    RotorOD/2 -- config.DuctMaxDia/2 (RotorOD/2 - Bridge) puts it right at
    the mechanical bridge, which the pocket geometry stays just inside of."""
    arc(radius, 0, config.SectorAngle, group=ROTOR_STEEL_GROUP)


def draw_rotor_inner_split(radius):
    """Full-sector arc at `radius`, just inside the deepest magnet-pocket
    point, separating the fine-meshed steel around the magnets (above) from
    the coarse bulk rotor steel (below). Crosses no pocket; it cuts only the
    RotorID -> duct (apbc2) span of each radial sector cut."""
    arc(radius, 0, config.SectorAngle, group=ROTOR_STEEL_GROUP)


def draw_rotor_skin_split(radius):
    """Full-sector arc at `radius`, between RotorEdgeRadius and RotorOD/2,
    separating the very dense outer rotor skin directly facing the airgap
    (above) from the rest of the dense edge band (below). Like the edge
    split, it never touches the RotorOD sliding-band arc itself."""
    arc(radius, 0, config.SectorAngle, group=ROTOR_STEEL_GROUP)


# ---------------------------------------
# Block label placement (one interior point per closed region) and
# material assignment.
# ---------------------------------------
def rotor_steel_label_point():
    """Interior point for the rotor lamination steel: sector's angular
    midpoint (clear of both magnet legs, which sit near angle 0 and
    SectorAngle), radius midway between rotor bore and rotor OD."""
    r = (config.RotorID / 2 + config.RotorOD / 2) / 2
    a = config.SectorAngle / 2
    return r * cosd(a), r * sind(a)


def stator_steel_label_point(angle_offset=0):
    """Interior point in the tooth back-iron at the given tooth-centerline
    angle, at a radius beyond the deepest slot-bottom fillet so it can't
    land inside a slot cavity."""
    c = config
    r_inner = c.StatorID / 2 + c.ShoeHeight + c.ShoeRadius + c.SlotHeight
    r = (r_inner + c.StatorOD / 2) / 2
    return r * cosd(angle_offset), r * sind(angle_offset)


def stator_slot_label_point(angle_offset=0):
    """Interior point at slot mid-body (mid-height between the shoe fillet
    and the slot bottom), at the given slot-centerline angle -- for placing
    a winding coil-side label. Same radius formula as reference.txt's
    0.5*(x5+x7)."""
    c = config
    r = c.StatorID / 2 + c.ShoeHeight + c.ShoeRadius + c.SlotHeight / 2
    return r * cosd(angle_offset), r * sind(angle_offset)


def setup_materials():
    """Import/define every material referenced by assign_block_labels(),
    using the exact names and property values from reference.txt (the
    working FEMM reference script for this motor)."""
    femm.mi_getmaterial("Air")

    femm.mi_addmaterial(
        "19 AWG", 1.0, 1.0, 0, 0, config.WireConductivity, 0, 0, 0, 3, 0, 0, 1,
        config.WireDiameter,
    )

    femm.mi_addmaterial("N36Z_20", 1.03, 1.03, 920000, 0, 0.667)

    femm.mi_addmaterial("M19_29G", 0, 0, 0, 0, 1.9, 0.34, 0, 0.94, 0)
    b_points = [
        0, 0.05, 0.1, 0.15, 0.36, 0.54, 0.65, 0.99, 1.2, 1.28, 1.33, 1.36,
        1.44, 1.52, 1.58, 1.63, 1.67, 1.8, 1.9, 2, 2.1, 2.3, 2.5,
        2.563994494, 3.779889874,
    ]
    h_points = [
        0, 22.28, 25.46, 31.83, 47.74, 63.66, 79.57, 159.15, 318.3, 477.46,
        636.61, 795.77, 1591.5, 3183, 4774.6, 6366.1, 7957.7, 15915, 31830,
        111407, 190984, 350135, 509252, 560177.2, 1527756,
    ]
    for b, h in zip(b_points, h_points):
        femm.mi_addbhpoint("M19_29G", b, h)


def assign_block_labels(split_stator=False, split_rotor=False):
    """Place one block label per closed region and assign its material.
    `split_stator`/`split_rotor` label the tooth/yoke and rotor-edge/bulk
    mesh splits (see draw_stator_tooth_yoke_split/draw_rotor_edge_split)
    instead of a single region -- the caller must have already drawn the
    corresponding split geometry."""
    c = config

    # NOTE: no per-slot "Air" label here. The slot mouths (SlotOpen, at the
    # stator bore) are left open -- no arc closes them off from the main
    # airgap annulus -- so all SectorSlots slot cavities and the airgap ring
    # are one single connected air region. A single label (the airgap one,
    # below) covers that whole region; adding one label per slot on top of
    # it puts multiple block labels in the same closed region, which FEMM
    # rejects ("more than one block label") at analysis time.

    if split_stator:
        # Below ToothSplitRadius, all SectorSlots+1 teeth (including the 2
        # half-teeth straddling the sector cuts) still touch each other at
        # the pinch points where neighboring slots' bottom fillets meet, so
        # they remain ONE connected region -- a single label (like the
        # unsplit case) at any interior tooth's centerline angle covers all
        # of them.
        r_tooth = (c.ShoeSplitRadius + c.ToothSplitRadius) / 2
        a_tooth = c.ToothPitch
        x, y = r_tooth * cosd(a_tooth), r_tooth * sind(a_tooth)
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop("M19_29G", 0, c.ToothMeshSize, "<None>", 0, STATOR_STEEL_GROUP, 0)
        femm.mi_clearselected()

        # Below ShoeSplitRadius, the slot mouths separate every tooth tip
        # into its own region: SectorSlots-1 full teeth at k*ToothPitch plus
        # the 2 half-teeth on the sector cuts (label nudged ToothPitch/8
        # inside the sector, still well within the tooth tip).
        r_shoe = (c.StatorID / 2 + c.ShoeSplitRadius) / 2
        shoe_angles = [c.ToothPitch / 8] + [
            k * c.ToothPitch for k in range(1, c.SectorSlots)
        ] + [c.SectorAngle - c.ToothPitch / 8]
        for a_shoe in shoe_angles:
            x, y = r_shoe * cosd(a_shoe), r_shoe * sind(a_shoe)
            femm.mi_addblocklabel(x, y)
            femm.mi_selectlabel(x, y)
            femm.mi_setblockprop("M19_29G", 0, c.ShoeMeshSize, "<None>", 0, STATOR_STEEL_GROUP, 0)
            femm.mi_clearselected()

        # The shoe arc also crosses every slot cavity, cutting a thin air
        # strip off the bottom of each slot body: bounded below by the
        # x3-x3m chord that closes the slot at shoe height (r ~3.232 in,
        # see build_model), above by ShoeSplitRadius, and on the sides by
        # the shoe/slot walls. One Air label per slot, at the slot center.
        x3 = c.StatorID / 2 + c.ShoeHeight
        y3 = (c.StatorPitchAirgap - c.SlotOpen) / 2
        r_chord = math.hypot(x3, y3) * cosd(c.ToothPitch / 2 - math.degrees(math.atan2(y3, x3)))
        r_strip = (r_chord + c.ShoeSplitRadius) / 2
        for k in range(c.SectorSlots):
            a_slot = c.StatorBaseRotation + k * c.ToothPitch
            x, y = r_strip * cosd(a_slot), r_strip * sind(a_slot)
            femm.mi_addblocklabel(x, y)
            femm.mi_selectlabel(x, y)
            femm.mi_setblockprop("Air", 0, c.ShoeMeshSize, "<None>", 0, SLOT_GROUP, 0)
            femm.mi_clearselected()

        r_yoke = (c.ToothSplitRadius + c.StatorOD / 2) / 2
        a_yoke = c.SectorAngle / 2
        x, y = r_yoke * cosd(a_yoke), r_yoke * sind(a_yoke)
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop("M19_29G", 0, c.YokeMeshSize, "<None>", 0, STATOR_STEEL_GROUP, 0)
        femm.mi_clearselected()
    else:
        x, y = stator_steel_label_point(c.StatorBaseRotation)
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop("M19_29G", 0, c.StatorSteelMeshSize, "<None>", 0, STATOR_STEEL_GROUP, 0)
        femm.mi_clearselected()

    if split_rotor:
        # Bulk: a clean annulus between RotorID and RotorInnerRadius, below
        # every magnet pocket, labelled at the sector's angular midpoint.
        a_bulk = c.SectorAngle / 2
        r_bulk = (c.RotorID / 2 + c.RotorInnerRadius) / 2
        x, y = r_bulk * cosd(a_bulk), r_bulk * sind(a_bulk)
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop("M19_29G", 0, c.RotorBulkMeshSize, "<None>", 0, ROTOR_STEEL_GROUP, 0)
        femm.mi_clearselected()

        # Steel around the magnets, RotorInnerRadius -> RotorEdgeRadius.
        # Both pockets lie entirely inside this layer; the steel above
        # each leg (pole shoe) stays connected to the rest through the thin
        # gap between the end-cap top (~3.096 in) and RotorEdgeRadius, so
        # the whole layer is still ONE region -- one label at the angular
        # midpoint (between the two legs) covers it.
        r_mid = (c.RotorInnerRadius + c.RotorEdgeRadius) / 2
        x, y = r_mid * cosd(a_bulk), r_mid * sind(a_bulk)
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop("M19_29G", 0, c.RotorMidMeshSize, "<None>", 0, ROTOR_STEEL_GROUP, 0)
        femm.mi_clearselected()

        # Dense edge band + very dense outer skin: two clean full-sector
        # annuli (RotorEdgeRadius -> RotorSkinRadius -> RotorOD), clear of
        # the magnet pockets, so one label per ring (angular midpoint)
        # covers each whole ring.
        a_edge = c.SectorAngle / 2
        for r0, r1, mesh_size in (
            (c.RotorEdgeRadius, c.RotorSkinRadius, c.RotorEdgeMeshSize),
            (c.RotorSkinRadius, c.RotorOD / 2, c.RotorSkinMeshSize),
        ):
            r_edge = (r0 + r1) / 2
            x, y = r_edge * cosd(a_edge), r_edge * sind(a_edge)
            femm.mi_addblocklabel(x, y)
            femm.mi_selectlabel(x, y)
            femm.mi_setblockprop("M19_29G", 0, mesh_size, "<None>", 0, ROTOR_STEEL_GROUP, 0)
            femm.mi_clearselected()
    else:
        x, y = rotor_steel_label_point()
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop("M19_29G", 0, c.RotorSteelMeshSize, "<None>", 0, ROTOR_STEEL_GROUP, 0)
        femm.mi_clearselected()

    # Airgap ring between the sliding band and the stator bore -- now closed
    # by the radial cuts added in draw_sector_outline. Manual mesh size is
    # required here: smartmesh's automesh=1 was verified to badly
    # undercompute torque (~53 vs the correct ~208 N*m for the same rotor
    # angle/currents) because it doesn't resolve the field across the
    # physically thin airgap finely enough.
    r = (SLIDING_BAND_RADIUS + c.StatorID / 2) / 2
    a = c.SectorAngle / 2
    x, y = r * cosd(a), r * sind(a)
    femm.mi_addblocklabel(x, y)
    femm.mi_selectlabel(x, y)
    femm.mi_setblockprop("Air", 0, c.AirgapMeshSize, "<None>", 0, AIRGAP_GROUP, 0)
    femm.mi_clearselected()

    for angle_offset, mirror in ((0, False), (c.SectorAngle, True)):
        duct_pt, magnet_pt, endcap_pt = rotor_pole_leg_label_points(angle_offset, mirror)
        magdir = rotor_pole_leg_magdirection(angle_offset, mirror)

        for pt, mesh_size in ((duct_pt, c.DuctMeshSize), (endcap_pt, c.MagnetAirMeshSize)):
            x, y = pt
            femm.mi_addblocklabel(x, y)
            femm.mi_selectlabel(x, y)
            femm.mi_setblockprop("Air", 0, mesh_size, "<None>", 0, MAGNET_GROUP, 0)
            femm.mi_clearselected()

        x, y = magnet_pt
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop("N36Z_20", 0, c.MagnetMeshSize, "<None>", magdir, MAGNET_GROUP, 0)
        femm.mi_clearselected()


def setup_circuits(time=0):
    """Add the 3-phase circuit properties (reference.txt's A/B/C), each
    carrying its instantaneous sinusoidal current at the given time -- 0 by
    default, matching reference.txt's static single-instant setup."""
    c = config
    for name, phase_shift in (("A", 0), ("B", 120), ("C", 240)):
        current = c.Current * sind(360 * c.Freq * time + c.Phase + phase_shift)
        femm.mi_addcircprop(name, current, 1)


def assign_windings():
    """One coil-side block label per stator slot in the sector (reference.txt's
    per-slot '19 AWG' winding assignment), using config.SlotCircuits/
    SlotCoilDirs -- one pole's worth of the machine's repeating AABBCC
    pattern; the sector's anti-periodic boundaries reproduce the opposite-
    direction half for the neighboring pole automatically.

    Slot centers sit at StatorBaseRotation + k*ToothPitch, k = 0..
    SectorSlots-1 -- the same angles tooth_centerline_point(0) resolves to
    for k=0 (ToothPitch/2 off the x-axis), since the half-unit + mirror +
    copyrotate construction in build_model tiles complete slots symmetrically
    within the sector, none touching the sector-cut boundaries."""
    c = config
    assert len(c.SlotCircuits) == c.SectorSlots == len(c.SlotCoilDirs), \
        "SlotCircuits/SlotCoilDirs must have exactly SectorSlots entries"

    for k in range(c.SectorSlots):
        angle_offset = c.StatorBaseRotation + k * c.ToothPitch
        x, y = stator_slot_label_point(angle_offset)
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        # automesh=1 (smartmesh-driven) for now -- manual c.WindingMeshSize
        # commented out: femm.mi_setblockprop("19 AWG", 0, c.WindingMeshSize, c.SlotCircuits[k], 0, SLOT_GROUP, c.SlotCoilDirs[k] * c.Turns)
        femm.mi_setblockprop(
            "19 AWG", 1, 0, c.SlotCircuits[k], 0, SLOT_GROUP,
            c.SlotCoilDirs[k] * c.Turns,
        )
        femm.mi_clearselected()


def assign_slot_air():
    """One Air block label per stator slot, coarsely meshed -- used in place
    of assign_windings() for analyses (e.g. cogging torque) that run
    unexcited, where the slot body carries no current and its exact mesh
    doesn't affect the field solution."""
    c = config
    for k in range(c.SectorSlots):
        angle_offset = c.StatorBaseRotation + k * c.ToothPitch
        x, y = stator_slot_label_point(angle_offset)
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop("Air", 0, c.SlotAirMeshSize, "<None>", 0, SLOT_GROUP, 0)
        femm.mi_clearselected()


# ---------------------------------------
# Boundary conditions: outer Dirichlet (A=0), anti-periodic sector cuts, and
# the anti-periodic sliding band across the airgap. The sector spans exactly
# one pole pitch (SectorAngle = 360/Npoles), so the sector's mirror image
# across each radial cut -- and across the sliding band -- carries opposite
# magnetic polarity, hence anti-periodic (not periodic) boundaries.
# ---------------------------------------
def setup_boundary_conditions(split_stator_radius=None, split_rotor_radius=None):
    """`split_stator_radius`/`split_rotor_radius` must match whatever was
    passed to draw_stator_tooth_yoke_split/draw_rotor_edge_split (or be
    left None if those weren't called) so the anti-periodic pairing below
    accounts for the radial segments those cuts auto-split."""
    c = config
    mid = c.SectorAngle / 2
    band_r = SLIDING_BAND_RADIUS

    # BdryFormat codes: 4 = Periodic, 5 = Anti-Periodic, 6 = Periodic air
    # gap, 7 = Anti-Periodic air gap. config.UsePeriodicBoundary switches
    # both the radial-cut pairs and the sliding band between the two, for
    # A/B comparison -- anti-periodic (the default) is correct when
    # SectorAngle spans an odd number of pole pitches, as it does here.
    radial_format = 4 if c.UsePeriodicBoundary else 5
    band_format = 6 if c.UsePeriodicBoundary else 7

    # Outer Dirichlet A=0 -- true exterior of the modeled domain.
    femm.mi_addboundprop("A=0", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    femm.mi_selectarcsegment(c.StatorOD / 2 * cosd(mid), c.StatorOD / 2 * sind(mid))
    femm.mi_setarcsegmentprop(1, "A=0", 0, 0)
    femm.mi_clearselected()

    # Anti-periodic pairs across the two radial sector cuts (angle 0 <->
    # SectorAngle). The RotorID -> RotorOD cut is not one continuous
    # segment: the rotor pole leg's duct nodes (x11, x12 in
    # _rotor_pole_leg_points, both on y=0/the sector cut line at
    # angle_offset=0) split it into 3 pieces where FEMM auto-splits the
    # segment at the intersecting nodes. So there are 5 segments total along
    # each sector cut line, each getting its own anti-periodic pairing:
    # stator body, rotor bore -> duct, duct itself, duct -> rotor OD, and
    # the band-to-stator airgap ring. The rotor-to-band side of the airgap
    # has no radial cut (left angularly open, part of the rotating mesh) --
    # its periodicity is handled by the sliding band's own anti-periodic
    # property (mySlidingBand) below.
    r11 = c.DuctMinDia / 2
    r12 = r11 + c.DuctThick / sind(c.alpha)
    # draw_stator_tooth_yoke_split's full-sector arc crosses both radial cut
    # lines, auto-splitting the single StatorID->StatorOD segment on each
    # side into 2 -- pair apbc1 as two ranges instead of one so each pairing
    # still matches exactly 2 segments (FEMM errors otherwise: "An
    # (anti)periodic BC is assigned to more than two segments").
    # draw_stator_shoe_split's arc adds one more such cut, so
    # `split_stator_radius` may be a single radius or a sequence of them;
    # the stator radial line is paired piecewise between consecutive cuts.
    if split_stator_radius is None:
        stator_radial_pairs = (("apbc1", c.StatorID / 2, c.StatorOD / 2),)
    else:
        if isinstance(split_stator_radius, (int, float)):
            split_stator_radius = (split_stator_radius,)
        stator_radii = [c.StatorID / 2] + sorted(split_stator_radius) + [c.StatorOD / 2]
        stator_radial_pairs = tuple(
            (f"apbc1{chr(ord('a') + i)}", r0, r1)
            for i, (r0, r1) in enumerate(zip(stator_radii[:-1], stator_radii[1:]))
        )
    # The rotor split arcs likewise cross both radial cut lines:
    # draw_rotor_edge_split/draw_rotor_skin_split land within the
    # r12->RotorOD/2 (apbc4) span, draw_rotor_inner_split within the
    # RotorID/2->r11 (apbc2) span. Each span is paired piecewise between
    # whichever cuts fall inside it (single radius or sequence).
    if split_rotor_radius is None:
        split_rotor_radius = ()
    elif isinstance(split_rotor_radius, (int, float)):
        split_rotor_radius = (split_rotor_radius,)

    def span_pairs(name, r_lo, r_hi):
        cuts = sorted(r for r in split_rotor_radius if r_lo < r < r_hi)
        if not cuts:
            return ((name, r_lo, r_hi),)
        radii = [r_lo] + cuts + [r_hi]
        return tuple(
            (f"{name}{chr(ord('a') + i)}", r0, r1)
            for i, (r0, r1) in enumerate(zip(radii[:-1], radii[1:]))
        )

    radial_pairs = stator_radial_pairs + span_pairs("apbc2", c.RotorID / 2, r11) + (
        ("apbc3", r11, r12),
    ) + span_pairs("apbc4", r12, c.RotorOD / 2) + (
        ("apbc5", band_r, c.StatorID / 2),
    )
    for name, r0, r1 in radial_pairs:
        femm.mi_addboundprop(name, 0, 0, 0, 0, 0, 0, 0, 0, radial_format, 0, 0)
        rmid = (r0 + r1) / 2
        for angle in (0, c.SectorAngle):
            femm.mi_selectsegment(rmid * cosd(angle), rmid * sind(angle))
        femm.mi_setsegmentprop(name, 0, 1, 0, 0)
        femm.mi_clearselected()

    # Sliding band at mid-airgap radius -- anti-periodic air-gap type
    # (BdryFormat 7 by default, matching antunes.fem's working
    # "mySlidingBand"; 6 for the periodic A/B comparison). Per FEMM's manual,
    # the last two params for an air-gap BdryFormat are "ia"/"oa" -- the
    # inner and outer boundary angles -- NOT an angle-plus-sector-span pair.
    # Both start at 0 (no rotation yet); only "ia" (propnum 10) is updated
    # via mi_modifyboundprop during a sweep. Previously "oa" was set to
    # c.SectorAngle here by mistake, permanently offsetting the outer
    # boundary's angle reference for the whole run.
    femm.mi_addboundprop(SLIDING_BAND_NAME, 0, 0, 0, 0, 0, 0, 0, 0, band_format, 0, 0)
    femm.mi_selectarcsegment(band_r * cosd(mid), band_r * sind(mid))
    femm.mi_setarcsegmentprop(1, SLIDING_BAND_NAME, 0, BAND_GROUP)
    femm.mi_clearselected()

    # RotorOD arc also carries the sliding band property, same group, so it
    # moves together with the mid-gap band when the rotor is rotated.
    # draw_rotor_edge_split's cut is a full-sector arc at a radius BELOW
    # RotorOD, so -- unlike the radial cuts -- it never touches the RotorOD
    # arc itself, which stays whole.
    femm.mi_selectarcsegment(c.RotorOD / 2 * cosd(mid), c.RotorOD / 2 * sind(mid))
    femm.mi_setarcsegmentprop(1, SLIDING_BAND_NAME, 0, BAND_GROUP)
    femm.mi_clearselected()


def pause(message):
    """Zoom to fit what's drawn so far and wait for the user to inspect it."""
    femm.mi_zoomnatural()
    try:
        input(f"{message} Press Enter to continue...")
    except EOFError:
        pass


def build_model(windings=True, split_stator_mesh=False, split_rotor_mesh=False):
    """Build the full sector geometry, materials, and boundary conditions
    into a fresh FEMM document. Leaves the document open and un-meshed
    (no mi_analyze) so callers can save/mesh/solve it as needed -- shared
    by the interactive __main__ below and by other scripts (e.g. a cogging
    torque sweep) that need the same model without the interactive pauses.

    windings=False fills the slots with plain Air (assign_slot_air) instead
    of the 19 AWG coil labels/circuits -- for unexcited analyses (cogging
    torque) where the copper isn't carrying current anyway.

    split_stator_mesh/split_rotor_mesh cut the stator steel into a dense
    tooth-shoe band + medium tooth body + coarse yoke, and the rotor steel into a dense edge band
    (facing the airgap) + coarse bulk, respectively -- see
    draw_stator_tooth_yoke_split/draw_rotor_edge_split.
    """
    femm.newdocument(0)

    femm.mi_probdef(0, "inches", "planar", 1e-8, config.Length, 10, 0)
    # smartmesh only affects blocks left on automesh=1; every block label in
    # assign_block_labels() below sets automesh=0 with its own numeric mesh
    # size (config.py), so each region (stator/rotor steel, airgap, magnet,
    # magnet-pocket air) meshes independently at its own resolution.
    femm.smartmesh(1)

    draw_sector_outline()

    draw_stator_half_unit(0)

    # Mirror-copy group 1 about the line through (0,0) and the tooth
    # centerline point -- ToothPitch/2 off the slot center -- to build the
    # neighboring half-tooth that would belong to the previous slot's
    # territory, leaving the original in place.
    xCutR, yCutR = tooth_centerline_point(0)
    femm.mi_selectgroup(SLOT_GROUP)
    femm.mi_mirror2(0, 0, xCutR, yCutR, 4)
    femm.mi_clearselected()

    # The half-unit draws only its own shoe-start point (x3); the mirror
    # produces the matching point on the neighboring half-tooth, but no
    # segment closes the tooth tip between them, leaving each slot cavity
    # open into the airgap. Add that closing segment once here, tagged
    # SLOT_GROUP so copyrotate2 below replicates it across every tooth.
    x3 = config.StatorID / 2 + config.ShoeHeight
    y3 = (config.StatorPitchAirgap - config.SlotOpen) / 2
    x3m, y3m = mirror_about_line(x3, y3, xCutR, yCutR)
    femm.mi_addsegment(x3, y3, x3m, y3m)
    femm.mi_selectsegment((x3 + x3m) / 2, (y3 + y3m) / 2)
    femm.mi_setgroup(SLOT_GROUP)
    femm.mi_clearselected()

    # The half-unit plus its mirror now form one complete slot spanning a
    # full ToothPitch. Copy-rotate that whole group about the origin,
    # ToothPitch per step, (SectorSlots - 1) more copies -- SectorSlots
    # slots total across the sector. Works for any Nslots/Npoles choice in
    # config.py, not just the current 48/8 -> 6-slot sector.
    femm.mi_selectgroup(SLOT_GROUP)
    femm.mi_copyrotate2(0, 0, config.ToothPitch, config.SectorSlots - 1, 4)
    femm.mi_clearselected()

    draw_rotor_pole_leg(0)

    draw_rotor_pole_leg(config.SectorAngle, mirror=True)

    if split_stator_mesh:
        draw_stator_shoe_split(config.ShoeSplitRadius)
        draw_stator_tooth_yoke_split(config.ToothSplitRadius)
    if split_rotor_mesh:
        draw_rotor_edge_split(config.RotorEdgeRadius)
        draw_rotor_skin_split(config.RotorSkinRadius)
        draw_rotor_inner_split(config.RotorInnerRadius)

    setup_materials()
    assign_block_labels(split_stator=split_stator_mesh, split_rotor=split_rotor_mesh)
    if windings:
        setup_circuits()
        assign_windings()
    else:
        assign_slot_air()
    setup_boundary_conditions(
        split_stator_radius=(
            (config.ShoeSplitRadius, config.ToothSplitRadius) if split_stator_mesh else None
        ),
        split_rotor_radius=(
            (config.RotorInnerRadius, config.RotorEdgeRadius, config.RotorSkinRadius)
            if split_rotor_mesh else None
        ),
    )


if __name__ == "__main__":
    femm.openfemm()
    build_model()

    pause("Geometry, materials, and boundary conditions complete.")

    femm.mi_saveas("ToyotaPrius.FEM")

    # Solve the saved model, then load the solution into the postprocessor
    # so mo_* calls (flux, torque, etc.) have something to act on.
    femm.mi_analyze()
    femm.mi_loadsolution()

    pause("Solution loaded. Manually mesh the model if needed.")
