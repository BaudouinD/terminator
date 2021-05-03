import math
from constants import *
from scipy.optimize import minimize
import numpy as np
import interpolation




# Given the sizes (a, b, c) of the 3 sides of a triangle, returns the angle between a and b using the alKashi theorem.
def alKashi(a, b, c, sign=-1):
    if a * b == 0:
        print("WARNING a or b is null in AlKashi")
        return 0
    # Note : to get the other altenative, simply change the sign of the return :
    return sign * math.acos(min(1, max(-1, (a ** 2 + b ** 2 - c ** 2) / (2 * a * b))))


def rotaton_2D(x, y, z, leg_angle):
    return rot_Z([x,y,z],leg_angle)


def rot_Y(pos, angle):
    return [np.cos(angle)*pos[0]-np.sin(angle)*pos[2],pos[1],np.sin(angle)*pos[0]+np.cos(angle)*pos[2]]

def rot_Z(pos, angle):
    return [np.cos(angle)*pos[0]-np.sin(angle)*pos[1],np.sin(angle)*pos[0]+np.cos(angle)*pos[1],pos[2]]

def rot_X(pos, angle):
    return [pos[0],np.cos(angle)*pos[1]-np.sin(angle)*pos[2],np.sin(angle)*pos[1]+np.cos(angle)*pos[2]]

def rad_to_deg(angle):
    return (angle * 180 )/ np.pi

def deg_to_rad(angle):
    return (angle * np.pi) / 180

# Computes the direct kinematics of a leg in the leg's frame
# Given the angles (theta1, theta2, theta3) of a limb with 3 rotational axes separated by the distances (l1, l2, l3),
# returns the destination point (x, y, z)
def computeDK(
    theta1,
    theta2,
    theta3,
    l1=constL1,
    l2=constL2,
    l3=constL3,
    use_rads=USE_RADS_INPUT,
    use_mm=USE_MM_OUTPUT,
):
    angle_unit = 1
    dist_unit = 1
    if not (use_rads):
        angle_unit = math.pi / 180.0
    if use_mm:
        dist_unit = 1000
    theta1 = THETA1_MOTOR_SIGN * theta1 * angle_unit
    theta2 = (THETA2_MOTOR_SIGN * theta2 - theta2Correction) * angle_unit
    theta3 = (THETA3_MOTOR_SIGN * theta3 - theta3Correction) * angle_unit

    planContribution = l1 + l2 * math.cos(theta2) + l3 * math.cos(theta2 + theta3)

    x = math.cos(theta1) * planContribution * dist_unit
    y = math.sin(theta1) * planContribution * dist_unit
    z = -(l2 * math.sin(theta2) + l3 * math.sin(theta2 + theta3)) * dist_unit

    return [x, y, z]


def computeDKDetailed(
    theta1,
    theta2,
    theta3,
    l1=constL1,
    l2=constL2,
    l3=constL3,
    use_rads=USE_RADS_INPUT,
    use_mm=USE_MM_OUTPUT,
):
    theta1_verif = theta1
    theta2_verif = theta2
    theta3_verif = theta3
    angle_unit = 1
    dist_unit = 1
    if not (use_rads):
        angle_unit = math.pi / 180.0
    if use_mm:
        dist_unit = 1000
    theta1 = THETA1_MOTOR_SIGN * theta1 * angle_unit
    theta2 = (THETA2_MOTOR_SIGN * theta2 - theta2Correction) * angle_unit
    theta3 = (THETA3_MOTOR_SIGN * theta3 - theta3Correction) * angle_unit

    planContribution = l1 + l2 * math.cos(theta2) + l3 * math.cos(theta2 + theta3)

    x = math.cos(theta1) * planContribution
    y = math.sin(theta1) * planContribution
    z = -(l2 * math.sin(theta2) + l3 * math.sin(theta2 + theta3))

    p0 = [0, 0, 0]
    p1 = [l1 * math.cos(theta1) * dist_unit, l1 * math.sin(theta1) * dist_unit, 0]
    p2 = [
        (l1 + l2 * math.cos(theta2)) * math.cos(theta1) * dist_unit,
        (l1 + l2 * math.cos(theta2)) * math.sin(theta1) * dist_unit,
        -l2 * math.sin(theta2) * dist_unit,
    ]
    p3 = [x * dist_unit, y * dist_unit, z * dist_unit]
    p3_verif = computeDK(
        theta1_verif, theta2_verif, theta3_verif, l1, l2, l3, use_rads, use_mm
    )
    if (p3[0] != p3_verif[0]) or (p3[1] != p3_verif[1]) or (p3[2] != p3_verif[2]):
        print(
            "ERROR: the DK function is broken!!! p3 = {}, p3_verif = {}".format(
                p3, p3_verif
            )
        )

    return [p0, p1, p2, p3]


# Computes the inverse kinematics of a leg in the leg's frame
# Given the destination point (x, y, z) of a limb with 3 rotational axes separated by the distances (l1, l2, l3),
# returns the angles to apply to the 3 axes
def computeIK(
    x,
    y,
    z,
    l1=constL1,
    l2=constL2,
    l3=constL3,
    verbose=False,
    use_rads=USE_RADS_OUTPUT,
    sign=-1,
    use_mm=USE_MM_INPUT,
):
    dist_unit = 1
    if use_mm:
        dist_unit = 0.001
    x = x * dist_unit
    y = y * dist_unit
    z = z * dist_unit

    # theta1 is simply the angle of the leg in the X/Y plane. We have the first angle we wanted.
    if y == 0 and x == 0:
        # Taking care of this singularity (leg right on top of the first rotational axis)
        theta1 = 0
    else:
        theta1 = math.atan2(y, x)

    # Distance between the second motor and the projection of the end of the leg on the X/Y plane
    xp = math.sqrt(x * x + y * y) - l1
    if xp < 0:
        print("Destination point too close")
        xp = 0

    # Distance between the second motor arm and the end of the leg
    d = math.sqrt(math.pow(xp, 2) + math.pow(z, 2))
    if d > l2 + l3:
        print("Destination point too far away")
        d = l2 + l3

    # Knowing l2, l3 and d, theta1 and theta2 can be computed using the Al Kashi law
    # There are 2 solutions for most of the points, forcing a convention here
    theta2 = alKashi(l2, d, l3, sign=sign) - Z_DIRECTION * math.atan2(z, xp)
    theta3 = math.pi + alKashi(l2, l3, d, sign=sign)

    if use_rads:

        result = [
            angleRestrict(THETA1_MOTOR_SIGN * theta1, use_rads=use_rads),
            angleRestrict(
                THETA2_MOTOR_SIGN * (theta2 + theta2Correction), use_rads=use_rads
            ),
            angleRestrict(
                THETA3_MOTOR_SIGN * (theta3 + theta3Correction), use_rads=use_rads
            ),
        ]

    else:
        result = [
            angleRestrict(THETA1_MOTOR_SIGN * math.degrees(theta1), use_rads=use_rads),
            angleRestrict(
                THETA2_MOTOR_SIGN * (math.degrees(theta2) + theta2Correction),
                use_rads=use_rads,
            ),
            angleRestrict(
                THETA3_MOTOR_SIGN * (math.degrees(theta3) + theta3Correction),
                use_rads=use_rads,
            ),
        ]
    if verbose:
        print(
            "Asked IK for x={}, y={}, z={}\n, --> theta1={}, theta2={}, theta3={}".format(
                x,
                y,
                z,
                result[0],
                result[1],
                result[2],
            )
        )

    return result

def angleRestrict(angle, use_rads=False):
    if use_rads:
        return modulopi(angle)
    else:
        return modulo180(angle)


# Takes an angle that's between 0 and 360 and returns an angle that is between -180 and 180
def modulo180(angle):
    if -180 < angle < 180:
        return angle

    angle = angle % 360
    if angle > 180:
        return -360 + angle

    return angle


def modulopi(angle):
    if -math.pi < angle < math.pi:
        return angle

    angle = angle % (math.pi * 2)
    if angle > math.pi:
        return -math.pi * 2 + angle

    return angle

 
def computeIKOriented(x,y,z,theta,leg_id,params,verbose=True):
    pos = (x, y, z*Z_DIRECTION)
    a = params.legAngles[leg_id-1]
    rot = np.array(rotaton_2D(x,y,z*Z_DIRECTION,a+theta))
    pos_ini = params.initLeg[leg_id-1] 
    if(verbose):
        print("rot :",rot)
        print("pos :", np.array(pos))
        print("pos apres rot :", np.array(pos) + np.array(pos_ini))
    res = []
    for i in range(3):
        res.append( rot[i] + np.array(pos_ini)[i])

    return computeIK(res[0],res[1],res[2])


def legs(leg1, leg2, leg3, leg4,leg5,leg6,theta,params):
    """w
    python simulator.py -m legs

    Le robot est figé en l'air, on contrôle toute les pattes

    - Sliders: les 12 coordonnées (x, y, z) du bout des 6 pattes
    - Entrée: des positions cibles (tuples (x, y, z)) pour le bout des 6 pattes
    - Sortie: un tableau contenant les 12 positions angulaires cibles (radian) pour les moteurs

    """
    ent = [leg1, leg2, leg3, leg4,leg5,leg6]
    targets = []
    for i in range(0,len(ent)):
        inv = computeIKOriented(ent[i][0], ent[i][1], ent[i][2],theta,i+1,params,verbose=False)
        for j in range(0,len(inv)):
            targets.append(inv[j])
    return targets

init_legs =  [[0.0,0.0,-0.06]] * 6
init_legs2 = [[0.0,0.0,0.0]] * 6
is_init = False
delta = 0
oldt = 0
splines = [interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D()]

def walk(t, speed_x, speed_y,theta,param,l1=constL1,
    l2=constL2,
    l3=constL3):
    """
    Le but est d'intégrer tout ce que nous avons vu ici pour faire marcher le robot
    - Sliders: speed_x, speed_y, speed_rotation, la vitesse cible du robot
    - Entrée: t, le temps (secondes écoulées depuis le début)
            speed_x, speed_y, et speed_rotation, vitesses cibles contrôlées par les sliders
    - Sortie: un tableau contenant les 12 positions angulaires cibles (radian) pour les moteurs
    """
    global init_legs2
    global is_init
    global delta
    global oldt
    global splines
    
    delta = t - oldt 
    if(is_init == False):
        print(init_legs2)
        is_init = True
    else:
        for i in range(0,len(splines)):
                splines[i].walk_trinalg(speed_x*delta, 
                                        speed_y*delta,
                                        param.z,l1,l2,l3,verbose=False)
                if(i%2):
                    init_legs2[i] = splines[i].interpolate(t%4)
                else:
                    init_legs2[i] = splines[i].interpolate((t+2)%4)
   
    oldt = t
    return legs(init_legs2[0],init_legs2[1],init_legs2[2],init_legs2[3],init_legs2[4],init_legs2[5],theta,param)

def legs_toupie(ent,params):
    targets = []
    for i in range(0,len(ent)):
        inv = computeIK(ent[i][0], ent[i][1], ent[i][2])
        for j in range(0,len(inv)):
            targets.append(inv[j])
    return targets

def toupie(t,speed,param,l1=constL1,l2=constL2,l3=constL3):
    global init_legs2
    global is_init
    global delta
    global oldt
    global splines
    
    delta = t - oldt 
    if(is_init == False):
        print(init_legs2)
        is_init = True
    else:
        for i in range(0,len(splines)):
                splines[i].toupie_trinalg(param.z,speed,l1,l2,l3,verbose=False)
                if(i%2):
                    init_legs2[i] = splines[i].interpolate(t%4)
                else:
                    init_legs2[i] = splines[i].interpolate((t+2)%4)

    oldt = t
    return legs_toupie(init_legs2,param)
