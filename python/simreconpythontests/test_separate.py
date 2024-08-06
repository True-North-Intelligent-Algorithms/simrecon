import sys
import os

from numpy import add

def addpaths():
    sys.path.insert(1, r'./python')

def test_makematrix():
    from simreconpython.separate import makematrix
    sepmatrix = makematrix(5,3)
    print(sepmatrix)
    #assert makematrix() == 'make matrix'


if __name__ == '__main__':
    print(os.getcwd())
    addpaths()
    test_makematrix()