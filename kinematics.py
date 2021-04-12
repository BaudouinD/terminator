import numpy as np
import constants
import utils
import interpolation

# THETA3_MOTOR_SIGN = -1
# THETA2_MOTOR_SIGN = 1
# THETA1_MOTOR_SIGN = 1
l1,l2,l3=0.045,0.065,0.087

def rot_Y(pos, angle):
    return [np.cos(angle)*pos[0]-np.sin(angle)*pos[2],pos[1],np.sin(angle)*pos[0]+np.cos(angle)*pos[2]]

def rot_Z(pos, angle):
    return [np.cos(angle)*pos[0]-np.sin(angle)*pos[1],np.sin(angle)*pos[0]+np.cos(angle)*pos[1],pos[2]]



def computeDK(a, b, c):
    a*=THETA1_MOTOR_SIGN
    b=THETA2_MOTOR_SIGN*b -theta2Correction
    c=THETA3_MOTOR_SIGN*b -theta3Correction
    c=-c
    MR3=[l3,0 ,0]
    MR3rot=rot_Y(MR3,c)

    MR2=[l2+MR3rot[0],MR3rot[1],MR3rot[2]]
    MR2rot=rot_Y(MR2,b)

    MR1=[l1+MR2rot[0],MR2rot[1],MR2rot[2]]
    MR1rot=rot_Z(MR1,a)


    return MR1rot


def al_kashi(adj1,adj2,opp):
    res=(adj1*adj1+adj2*adj2-opp*opp)/(2*adj1*adj2)
    newres=min(max(res,-1),1)
    if(newres!=res):
       print("WARNING : value in arccos out of bounds ([-1;1]) :%f " %res)
    return np.arccos(newres)

def rotaton_2D(x, y, z, leg_angle):
    return rot_Z([x,y,z],leg_angle)

def computeIK(x,y,z):
    theta1=np.arctan2(y,x)
    AH=np.sqrt(x*x+y*y)-l1
    AM=np.sqrt(z*z+AH*AH)
    gamma2=al_kashi(l2,l3,AM)#np.arccos((l2*l2+l3*l3-AM*AM)/(2*l2*l3))
    gamma=np.pi-gamma2
    delta=al_kashi(l2,AM,l3)#np.arccos((l2*l2+AM*AM-l3*l3)/(2*l2*AM))
 
    alpha=np.arctan2(-z,AH)
    beta=delta-alpha

    return [theta1, beta, gamma]



def legs(leg1, leg2, leg3, leg4):
    """w
    python simulator.py -m legs

    Le robot est figé en l'air, on contrôle toute les pattes

    - Sliders: les 12 coordonnées (x, y, z) du bout des 4 pattes
    - Entrée: des positions cibles (tuples (x, y, z)) pour le bout des 4 pattes
    - Sortie: un tableau contenant les 12 positions angulaires cibles (radian) pour les moteurs

    """
    ent = [leg1, leg2, leg3, leg4]
    angle = [np.pi*5.0/4.0,np.pi*3.0/4.0,np.pi/4.0,-np.pi/4.0]
    targets = []
    for i in range(0,len(ent)):
        
        tmp = np.array(ent[i])
        tmp = MatRotation(angle[i],2).dot(tmp)
        tmp[0] = tmp[0] - 0.040  
        inv = computeIK(tmp[0], tmp[1], tmp[2])
        for j in range(0,len(inv)):
            targets.append(inv[j])
    return targets


#splines = [interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D()]i
init_legs = [[-0.1,0.1,-0.05],[-0.1,-0.1,-0.05],[0.1,-0.1,-0.05],[0.1,0.1,-0.05]]
init_legs2 = [[-0.1,0.1,-0.05],[-0.1,-0.1,-0.05],[0.1,-0.1,-0.05],[0.1,0.1,-0.05]]

def walk(t, speed_x, speed_y, speed_rotation):
    """
    python simulator.py -m walk

    Le but est d'intégrer tout ce que nous avons vu ici pour faire marcher le robot

    - Sliders: speed_x, speed_y, speed_rotation, la vitesse cible du robot
    - Entrée: t, le temps (secondes écoulées depuis le début)
            speed_x, speed_y, et speed_rotation, vitesses cibles contrôlées par les sliders
    - Sortie: un tableau contenant les 12 positions angulaires cibles (radian) pour les moteurs
    """
    global init_legs2
    targets = [0]*12


   
    splines = [interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D(),interpolation.LinearSpline3D()]
    for i in range(0,len(splines)):
        print("i=",i)
        splines[i].walk_trinalg((init_legs[i]),[speed_x+init_legs[i][0],speed_y+ init_legs[i][1],-0.05])
        if( i == 0 or i == 2):
            init_legs2[i] = splines[i].interpolate(t%4)
        else:
            init_legs2[i] = splines[i].interpolate((t+2)%4)

        

    return legs(init_legs2[0],init_legs2[1],init_legs2[2],init_legs2[3])

