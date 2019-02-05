import cmath
import math
import numpy as np
import scipy as sp

posi = complex(0, 1)
negi = complex(0, -1)

#Zero State |0>
zero = np.matrix([[1],
                [0]])

#One State |1>
one = np.matrix([[0],
                [1]])

#Plus State |+>
plus = 1/math.sqrt(2) * np.matrix([[1],
                                [1]])

#Minus State |->
minus = 1/math.sqrt(2) * np.matrix([[1],
                                    [-1]])

class gates:

    """ONE QBIT"""
    #NOT Gate
    def xgate(qinput):
        xmat = np.matrix([[0, 1],
                        [1, 0]])
        return np.dot(xmat, qinput)

    #Sqaure Root NOT Gate
    def sqrtxgate(qinput):
        sxmat = np.matrix([[posi, negi],
                        [negi, posi]])
        return np.dot(sxmat, qinput)

    #Z Gate
    def zgate(qinput):
        zmat = np.matrix([[1, 0],
                        [0, -1]])
        return np.dot(zmat, qinput)

    #Y Gate
    def ygate(qinput):
        ymat = np.matrix([[0, negi],
                        [posi, 0]])
        return np.dot(ymat, qinput)

    #Hadamard Gate
    def hadamard(qinput):
        hmat = 1/math.sqrt(2) * np.matrix([[1, 1],
                                        [1, -1]])
        return np.dot(hmat, qinput)

    """TWO QBITS"""
    #CNOT Gate
    def cnot(qbitone, qinput):
        cmat = np.matrix([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [1, 0, 0, 1],
                        [0, 1, 1, 0]])
        retarr = np.dot(cmat, qinput)
        np.put(retarr, [0, 1], [qbitone[0], qbitone[1]])
        return retarr

    #Swap Gate
    def swap(qinput):
        smat = np.matrix([[1, 0, 0, 0],
                        [0, 0, 1, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 1]])
        return np.dot(smat, qinput)

    #Square Root Swap
    def sqrtswap(qinput):
        ssmat = np.matrix([[1, 0, 0, 0],
                        [0 , 1/2 * posi, 1/2 * negi, 0],
                        [0, 1/2 * negi, 1/2 * posi, 0],
                        [0, 0, 0, 1]])
        return np.dot(ssmat, qinput)

    #Ising (XX) Gate
    def ising(qinput, qphi):
        imat = 1/math.sqrt(2) * np.matrix([[1, 0, 0, cmath.exp(posi*(1 - math.pi/2))],
                                        [0, 1, negi, 0],
                                        [0, negi, 1, 0],
                                        [cmath.exp(posi*(1 - math.pi/2)), 0, 0, 1]])
        return np.dot(imat, qinput)

    """THREE QBITS"""
    #Toffoli Gate
    def toffoli(qinput):
        tmat = np.matrix([[1, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 1, 0, 0, 0, 0, 0],
                        [0, 0, 0, 1, 0, 0, 0, 0],
                        [0, 0, 0, 0, 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 1, 0, 0],
                        [1, 0, 0, 0, 0, 0, 0, 1],
                        [1, 0, 0, 0, 0, 0, 1, 0]])
        return np.dot(tmat, qinput)

    #Fredkin Gate
    def fredkin(qinput):
        fmat = np.matrix([[1, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 1, 0, 0, 0, 0, 0],
                        [0, 0, 0, 1, 0, 0, 0, 0],
                        [0, 0, 0, 0, 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 1, 0],
                        [1, 0, 0, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 1]])
        return np.dot(fmat, qinput)

    #Peres Gate
    def peres(qinput):
        pmat = np.matrix([[1, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 1, 0, 0, 0, 0, 0],
                        [0, 0, 0, 1, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 0, 0, 0, 1],
                        [0, 0, 0, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 1, 0, 0, 0]])
        return np.dot(pmat, qinput)

    """INF QBITS"""
    #CU Gate
    def cu(qinput):
        qtemp = qinput
        print(qinput[[0]])
        if qtemp[[0]] == 1:
            for i in range(len(qtemp)):
                qtemp[i + 1] = gates.xgate(qinput[i + 1])
            return qtemp

class assembly:

    def nkron(*args):
        result = np.array([[1]])
        for op in args:
            result = np.kron(result, op)
        return result
print(one)
print()
print(zero)
print()
print(gates.swap(np.kron(zero, one)))
print(gates.swap(np.kron(one, zero)))
