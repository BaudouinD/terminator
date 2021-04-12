import numpy as np
import constants
import utils

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