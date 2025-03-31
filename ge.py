import argparse
import json
from sys import argv
from sage.all import *

class GePy:
    """
    main functionality, to perform all computations in Sage
    """
    def __init__(self, n, frame_matrix, complex_structure, connection="lc"):
        # setup init self variables
        self.n = n
        # setup manifold with coordinates and as parallelizable
        M = Manifold(Integer(self.n), 'M', field="real")
        x = tuple(f"x{i}" for i in range(self.n))
        U = M.chart(names=x)
        x = U._first_ngens(self.n) # makes M parallelizable
        # parse the orthonormal frame 
        frame = M.automorphism_field()
        frame[:] = matrix(self.n, self.n, lambda i, j: sage_eval(frame_matrix[i][j], locals={'x':x}))
        # orthonormal frame and coframe to self
        self.ed = M.default_frame().new_frame(frame, 'e') # orthonormal frame
        self.eu = self.ed.dual_basis() # orthonormal coframe 
        self.xd = M.frames()[0] # coordinate frame
        self.xu = self.xd.dual_basis() # ooordinate coframe
        # if complex manifold, parse complex structure, check integrability, compute fundamental 2-form, 
        if complex_structure:
            self.complex_parsing()


    def complex_parsing(self):
        # parse the complex structure
        self.J = M.tensor_field(1, 1, name='J')
        for i in range(self.n):
            for j in range(self.n):
                self.J[i, j] = complex_structure[i][j]
        print("Complex structure integrable:", self.is_integrable(self.J))
        # compute fundamental 2-form
        self.F = M.diff_form(Integer(2), name="F")
        for i in range(DIM):
            for j in range(DIM):
                F[ed, i, j] = g(J(ed[i]), ed[j])
        self.dF = F.derivative()


    def nijenhuis(self, X, Y):
        """
        return True if Nijenhuis tensor is 0, and False otherwise
        """
        return X.bracket(Y) + self.J(self.J(X).bracket(Y) + X.bracket(self.J(Y))) - self.J(X).bracket(Y)


    def is_integrable(self, J):
        """
        check if provided almost complex structure is integrable using Nijenhuis tensor
        """
        return all([self.nijenhuis(self.xd[i],self.xd[j]).contract(self.xu[k]) == 0 for i in range(self.n) for j in range(self.n) for k in range(self.n)])


def parse_arguments():
    """
    arguments parser using library argparse
    """
    parser = argparse.ArgumentParser(prog='gepy', description="Geometric computations in Python")

    parser.add_argument("jsonfile")
    parser.add_argument("-c", "--connection", choices=["levicivita","gauduchon"], type=str, help="Connection to compute and utilise, default is Levi-Civita")
    parser.add_argument("-t", type=int, help="Parameter for the Gauduchon connection, e.g. '-t 1' is Chern connection)")

    args = parser.parse_args()

    if args.connection == "gauduchon" and args.t is None:
        raise argparse.ArgumentError(None, "The '-t' argument is required when '-c gauduchon'.")
    if args.connection is None:
        args.connection = "levicivita"

    return args


def square_matrix(matrix):
    """
    returns dimension of matrix if square, 0 otherwise
    """
    n = len(matrix)
    for r in matrix:
        if len(r) != n:
            return 0
    return n


def parse_json(filename):
    """
    parse the json file containing orthonormal frame
    """
    with open(filename, "r") as f:
        data = json.load(f)
    # check if orthonormal frame in json, required
    try:
        frame = data["frame"]
    except KeyError:
        raise Exception("Orthonormal frame ('frame') required in json file")
    # check if complex structure in json, optional
    try:
        complex_structure = data["complex_structure"]
    except KeyError:
        complex_structure = None
    # the frame needs to be a square matrix, check len rows == len columns
    frame_dim = square_matrix(frame)
    complex_structure_dim = square_matrix(complex_structure)
    if frame_dim == 0 or complex_structure_dim == 0 or frame_dim != complex_structure_dim:
        raise Exception("Matrices in json must be square and have the same dimension")

    return frame_dim, frame, complex_structure


if __name__ == "__main__":
    ARGS = parse_arguments()
    DIM, FRAME, COMPLEX_STRUCTURE = parse_json(ARGS.jsonfile)
    GEPY = GePy(DIM, FRAME, COMPLEX_STRUCTURE, ARGS.connection)
    # TODO: implement all Ricci curvatures for Chern connection

