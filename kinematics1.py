import numpy as np
from constants import *
import utils
import interpolation

l1,l2,l3=constL3,constL2,constL1


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

def computeDK(a, b, c, use_rads=True):

    if(use_rads):
        a*=THETA1_MOTOR_SIGN
        b=THETA2_MOTOR_SIGN*b  - theta2Correction
        c=THETA3_MOTOR_SIGN*c  - theta3Correction
    else:
        a*=THETA1_MOTOR_SIGN
        b=THETA2_MOTOR_SIGN*b  - rad_to_deg(theta2Correction)
        c=THETA3_MOTOR_SIGN*c  - rad_to_deg(theta3Correction)

    MR3=[l3,0 ,0]
    MR3rot=rot_Y(MR3,c)

    MR2=[l2+MR3rot[0],MR3rot[1],MR3rot[2]]
    MR2rot=rot_Y(MR2,b)

    MR1=[l1+MR2rot[0],MR2rot[1],MR2rot[2]]
    MR1rot=rot_Z(MR1,a)

    return  list(map(lambda x : x* 1000,MR1rot)) #list(map(rad_to_deg,MR1rot))

print("DK ",computeDK(0,0,0))

def computeDKDetailed(a, b, c , use_rads=True):

    if(use_rads):
        a*=THETA1_MOTOR_SIGN
        b=THETA2_MOTOR_SIGN*b  - theta2Correction
        c=THETA3_MOTOR_SIGN*c  - theta3Correction
    else:
        a*=THETA1_MOTOR_SIGN
        b=THETA2_MOTOR_SIGN*b  - rad_to_deg(theta2Correction)
        c=THETA3_MOTOR_SIGN*c  - rad_to_deg(theta3Correction)

    MR3=[l3,0 ,0]
    MR3rot=rot_Y(MR3,c)

    MR2=[l2+MR3rot[0],MR3rot[1],MR3rot[2]]
    MR2rot=rot_Y(MR2,b)

    MR1=[l1+MR2rot[0],MR2rot[1],MR2rot[2]]
    MR1rot=rot_Z(MR1,a)
    
    return [[0,0,0],list(map(lambda x : x* 1000,MR3rot)),list(map(lambda x : x* 1000,MR2rot)),list(map(lambda x : x* 1000,MR1rot))]#[MR3rot,MR2rot,MR1rot,[0,0,0]]

print("Detailed ",computeDKDetailed(0,0,0))

def al_kashi(adj1,adj2,opp):
    res=(adj1*adj1+adj2*adj2-opp*opp)/(2*adj1*adj2)
    newres=min(max(res,-1),1)
    if(newres!=res):
       print("WARNING : value in arccos out of bounds ([-1;1]) :%f " %res)
    return np.arccos(newres)

def rotaton_2D(x, y, z, leg_angle):
    return rot_Z([x,y,z],leg_angle)


def computeIK(x,y,z,use_rads=True):

    z = -z

    theta1=np.arctan2(y,x)
    AH=np.sqrt(x*x+y*y)-l1
    AM=np.sqrt(z*z+AH*AH)
    gamma2=al_kashi(l2,l3,AM) #np.arccos((l2*l2+l3*l3-AM*AM)/(2*l2*l3))
    gamma=np.pi-gamma2 + theta3Correction
    delta=al_kashi(l2,AM,l3) #np.arccos((l2*l2+AM*AM-l3*l3)/(2*l2*AM))
 
    alpha=np.arctan2(-z,AH)
    beta=delta-alpha + theta2Correction

    if(use_rads):
        return [theta1, beta, gamma] 
    else:
        return list(map(rad_to_deg,[theta1, beta, gamma]))



def legs(leg1, leg2, leg3, leg4,leg5,leg6):
    """w
    python simulator.py -m legs

    Le robot est figé en l'air, on contrôle toute les pattes

    - Sliders: les 12 coordonnées (x, y, z) du bout des 6 pattes
    - Entrée: des positions cibles (tuples (x, y, z)) pour le bout des 6 pattes
    - Sortie: un tableau contenant les 12 positions angulaires cibles (radian) pour les moteurs

    """
    ent = [leg1, leg2, leg3, leg4,leg5,leg6]
    angle = [math.pi / 2, 0, 0, -math.pi / 2, math.pi, math.pi]
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
        if (i== 100):
            print("i=",i)
            splines[i].walk_trinalg(init_legs[i],[speed_x+init_legs[i][0],speed_y+ init_legs[i][1],-0.05])
            if(i%2):
                init_legs[i] = splines[i].interpolate(t%4)
            else:
                init_legs[i] = splines[i].interpolate((t+2)%4)

        

    return legs(init_legs[0],init_legs[1],init_legs[2],init_legs[3],init_legs[4],init_legs[5])

