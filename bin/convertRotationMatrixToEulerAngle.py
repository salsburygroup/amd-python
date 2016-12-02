#Convert Rotation Matrix to Euler Angle
#Ryan Godwin
#Modified from:https://www.learnopencv.com/rotation-matrix-to-euler-angles/

#Takes a rotation matrix from VMD output and converts to Euler Angles
#For example
#python /Users/fwamps/SalsburyGroup/Python/sandbox/convertRotationMatrixToEulerAngle.py
#-rMat '{{0.31456 0.948179 0.0448005} {0.785671 -0.286553 0.548277} {0.532702 -0.137269 -0.835096}}'

import numpy as np
import math
import re
import argparse

#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Convert a rotation matrix output from VMD into Euler Angles', add_help=False) 

#List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-rMat', action='store', dest='Rin',help='The Rotation Matrix from which to determine Euler Angles',type=str,required=True)

#Parse into useful form
UserInput=parser.parse_args()

#np.set_printoptions(precision=5)

Rin2=re.sub('[{}]','',UserInput.Rin)
RinClean=Rin2.split(' ')
RinStr=np.array(RinClean, dtype='f')
RinNum=RinStr.astype(np.float)
R=np.reshape(RinNum, (3,3))

print(R)

# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6
 
 
# Calculates rotation matrix to euler angles
# The result is the same as MATLAB except the order
# of the euler angles ( x and z are swapped ).
def rotationMatrixToEulerAngles(R) :
 
    #assert(isRotationMatrix(R))
     
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
     
    singular = sy < 1e-6
 
    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])
    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0
    return np.array([x*(180/math.pi), y*(180/math.pi), z*(180/math.pi)])

print(rotationMatrixToEulerAngles(R))