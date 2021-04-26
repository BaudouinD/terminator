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
    # print(
    #     "corrected angles={}, {}, {}".format(
    #         theta1 * (1.0 / angle_unit),
    #         theta2 * (1.0 / angle_unit),
    #         theta3 * (1.0 / angle_unit),
    #     )
    # )

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

    # print(
    #     "corrected angles={}, {}, {}".format(
    #         theta1 * (1.0 / angle_unit),
    #         theta2 * (1.0 / angle_unit),
    #         theta3 * (1.0 / angle_unit),
    #     )
    # )

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
    # if xp < 0:
    #     print("Destination point too close")
    #     xp = 0

    # Distance between the second motor arm and the end of the leg
    d = math.sqrt(math.pow(xp, 2) + math.pow(z, 2))
    # if d > l2 + l3:
    #     print("Destination point too far away")
    #     d = l2 + l3

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


# def al_kashi(adj1,adj2,opp):
#     res=(adj1*adj1+adj2*adj2-opp*opp)/(2*adj1*adj2)
#     newres=min(max(res,-1),1)
#     if(newres!=res):
#        print("WARNING : value in arccos out of bounds ([-1;1]) :%f " %res)
#     return np.arccos(newres)

# def computeIK(x,y,z, l1=constL1,l2=constL2,l3=constL3,use_rads=True):

#     z = -z

#     theta1=np.arctan2(y,x)
#     AH=np.sqrt(x*x+y*y)-l1
#     AM=np.sqrt(z*z+AH*AH)
#     gamma2=al_kashi(l2,l3,AM) #np.arccos((l2*l2+l3*l3-AM*AM)/(2*l2*l3))
#     gamma=np.pi-gamma2 + theta3Correction
#     delta=al_kashi(l2,AM,l3) #np.arccos((l2*l2+AM*AM-l3*l3)/(2*l2*AM))
 
#     alpha=np.arctan2(-z,AH)
#     beta=delta-alpha + theta2Correction

#     if(use_rads):
#         return [theta1, beta, gamma] 
#     else:
#         return list(map(rad_to_deg,[theta1, beta, gamma]))


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


def computeIKOriented(x,y,z,leg_id,params,verbose=True):
    print(x,y,z)
    x += params.initLeg[leg_id - 1][0]
    y += params.initLeg[leg_id - 1][1]
    z += params.z
   
    tab= rotaton_2D(x,y,z,params.legAngles[leg_id - 1])
    res = computeIK(tab[0],tab[1],tab[2])
    print(res)

    return res


def legs(leg1, leg2, leg3, leg4,leg5,leg6):
    """w
    python simulator.py -m legs

    Le robot est figé en l'air, on contrôle toute les pattes

    - Sliders: les 12 coordonnées (x, y, z) du bout des 6 pattes
    - Entrée: des positions cibles (tuples (x, y, z)) pour le bout des 6 pattes
    - Sortie: un tableau contenant les 12 positions angulaires cibles (radian) pour les moteurs

    """
    ent = [leg1, leg2, leg3, leg4,leg5,leg6]
    angle = [math.pi / 4,
        -math.pi / 4,
        -math.pi / 2,
        -3 * math.pi / 4,
        3 * math.pi / 4,
        math.pi / 2,]
    targets = []
    for i in range(0,len(ent)):
        
        tmp = np.array(ent[i])
        #tmp = MatRotation(angle[i],2).dot(tmp)
        tmp = rot_Z(tmp,angle[i])
        #tmp[0] = tmp[0] - 0.040  
        inv = computeIK(tmp[0], tmp[1], tmp[2])
        for j in range(0,len(inv)):
            targets.append(inv[j])
    return targets


#splines = [interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D()]i
init_legs = [[-0.1,0.1,-0.05],[-0.1,-0.1,-0.05],[0.1,-0.1,-0.05],[0.1,0.1,-0.05],[0.1,0.1,-0.05],[0.1,0.1,-0.05]]
init_legs2 = [[-0.1,0.1,-0.05],[-0.1,-0.1,-0.05],[0.1,-0.1,-0.05],[0.1,0.1,-0.05],[0.1,0.1,-0.05],[0.1,0.1,-0.05]]

def walk(t, speed_x, speed_y,param):
    """
    Le but est d'intégrer tout ce que nous avons vu ici pour faire marcher le robot
    - Sliders: speed_x, speed_y, speed_rotation, la vitesse cible du robot
    - Entrée: t, le temps (secondes écoulées depuis le début)
            speed_x, speed_y, et speed_rotation, vitesses cibles contrôlées par les sliders
    - Sortie: un tableau contenant les 12 positions angulaires cibles (radian) pour les moteurs
    """
    #global init_legs2
    #targets = [0]*12
    init_legs=param.initLeg
    print(init_legs)
    splines = [interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D()]
    for i in range(0,len(splines)):
        if (i==100):
            print("i=",i)
            splines[i].walk_trinalg(init_legs[i],[speed_x+init_legs[i][0],speed_y+ init_legs[i][1],-0.05])
            if(i%2):
                init_legs[i] = splines[i].interpolate(t%4)
            else:
                init_legs[i] = splines[i].interpolate((t+2)%4)

    

    return legs(init_legs[0],init_legs[1],init_legs[2],init_legs[3],init_legs[4],init_legs[5])
