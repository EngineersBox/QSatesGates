import numpy as np
import numpy.random
import math
import cmath
import scipy as sp
import scipy.linalg

posi = complex(0, 1)
negi = complex(0, -1)

normalisesate = lambda state: state / sp.linalg.norm(state)

#Equivalent of |1>
one = np.array([[0.0],
                  [1.0]])

#Equivalent of |0>
zero = np.array([[1.0],
                   [0.0]])

plus = normalisesate(zero + one)

class gates:

    def spos(qinputone, qinputtwo):
        return normalisesate(qinputone + qinputtwo)

    def xgate(qinput):
        xmat = np.array([[0, 1],
                          [1, 0]])
        return np.dot(xmat, qinput)

    def ygate(qinput):
        ymat = np.array([[0, negi],
                         [posi, 0]])
        return np.dot(ymat, qinput)

    def zgate(qinput):
        zmat = np.array([[1, 0],
                         [0, -1]])
        return np.dot(zmat, qinput)

    def hadamard(qinput):
        hmat = 1. / np.sqrt(2) * np.array([[1, 1],
                                           [1, -1]])
        return np.dot(hmat, qinput)

    def phase(qinput):
        pmat = np.array([[1, 0],
                         [0, posi]])
        return np.dot(pmat, qinput)

    def phaserot(qinput, theta):
        rmat = np.array([[1, 0],
                         [0, cmath.exp(posi * theta)]])
        return np.dot(rmat, qinput)

    def cnot(qinputone, qinputtwo):
        if qinputone[0] == 1 and qinputone[1] == 0:
            qtwo = gates.xgate(qinputtwo)
            return qtwo
        else:
            return qinputtwo

    def swap(qinputone, qinputtwo): #add for loop iteration here to swap elements
        qtemp = qinputone
        outone = qinputtwo
        outtwo = qtemp
        return outone, outtwo

    def nkron(*args):
        result = np.array([[1.0]])
        for op in args:
            result = np.kron(result, op)
        return result

class interaction:

    def measurement(qinputone, qinputtwo):
        catstate = normalisesate(qinputone + qinputtwo)
        rhocatstate = np.dot(catstate, catstate.T)

        p0 = np.dot(zero, zero.T)
        p1 = np.dot(one, one.T)
        Id = np.eye(2)
        prob0 = np.trace(np.dot(gates.nkron(p0, Id), rhocatstate))

        if np.random.rand() < prob0:
            result = 0
            resultstate = normalisesate(np.dot(gates.nkron(p0, Id), catstate))
        else:
            result = 1
            resultstate = normalisesate(np.dot(gates.nkron(p1, Id), catstate))

#qbit = np.dot(xgate, qbit)
print(gates.swap(one, zero))
