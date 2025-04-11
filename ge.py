import argparse
import json
import sage.all as sg

class GePy:
    """
    main functionality, to perform all computations in Sage
    """
    def __init__(self, n, frame_matrix):
        # setup init self variables
        self.n = n
        # setup manifold with coordinates and as parallelizable
        self.M = sg.Manifold(sg.Integer(self.n), "M", field="real")
        x = tuple(f"x{i}" for i in range(self.n))
        U = self.M.chart(names=x)
        x = U._first_ngens(self.n) # makes M parallelizable
        # parse the orthonormal frame 
        frame = self.M.automorphism_field()
        frame[:] = sg.matrix(self.n, self.n, lambda i, j: sg.sage_eval(frame_matrix[i][j], locals={'x':x}))
        # orthonormal frame and coframe to self
        self.ed = self.M.default_frame().new_frame(frame, 'e') # orthonormal frame
        self.eu = self.ed.dual_basis() # orthonormal coframe 
        self.xd = self.M.frames()[0] # coordinate frame
        self.xu = self.xd.dual_basis() # ooordinate coframe


    def get_connection(self, connection_name):
        assert type(self.nab) is sg.sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection
        if connection_name == "levicivita":
            return self.nab
        assert type(self.nab_t) is sg.sage.manifolds.differentiable.affine_connection.AffineConnection
        if connection_name == "gauduchon":
            return self.nab_t
        return None


    def compute_curvature(self, nabla):
        assert type(nabla) is sg.sage.manifolds.differentiable.affine_connection.AffineConnection or\
                sg.sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection
        return nabla.riemann()


    def compute_ricci_curvature(self, nabla):
        assert type(nabla) is sg.sage.manifolds.differentiable.affine_connection.AffineConnection or\
                sg.sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection
        return nabla.ricci()


    def compute_scalar_curvature(self, nabla):
        assert type(nabla) is sg.sage.manifolds.differentiable.affine_connection.AffineConnection or\
                sg.sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection
        return sum( self.g.inverse()[i,j]*nabla.ricci()[i,j] for i in range(self.n) for j in range(self.n) )


    def compute_sectional_curvature(self, nabla):
        assert type(nabla) is sg.sage.manifolds.differentiable.affine_connection.AffineConnection or\
                sg.sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection
        sectional_curvature = dict()
        for i in range(self.n):
            for j in range(i+1, self.n):
                sectional_curvature[(i,j)] = nabla.riemann()[i,j,j,i] / (self.g[i,i]*self.g[j,j] - self.g[i,j]**2)
        return sectional_curvature


    def compute_torsion(self, nabla):
        assert type(nabla) is sg.sage.manifolds.differentiable.affine_connection.AffineConnection or\
                sg.sage.manifolds.differentiable.levi_civita_connection.LeviCivitaConnection
        return nabla.torsion()


    def compute_metric(self):
        self.g = self.M.metric('g')
        self.g[self.ed,:] = sg.identity_matrix(self.n)


    def gauduchon_connection(self, t):
        self.nab_t = self.M.affine_connection("nab_t")

        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    self.nab_t[self.ed,i,j,k] = self.nab[self.ed,i,j,k] - ((1-t)/4)*self.dF(self.J(self.ed[i]),self.J(self.ed[j]),self.J(self.ed[k])) + \
                            ((1+t)/4)*self.dF(self.J(self.ed[i]),self.ed[j],self.ed[k]) - (1/2)*self.nijenhuis(self.xd[i],self.xd[j]).contract(self.xu[k])


    def levi_civita_connection(self):
        self.nab = self.g.connection()


    def complex_parsing(self, complex_structure):
        # parse the complex structure
        self.J = self.M.tensor_field(1, 1, name='J')
        for i in range(self.n):
            for j in range(self.n):
                self.J[i, j] = complex_structure[i][j]
        print("Complex structure integrable:", self.is_integrable())
        # compute fundamental 2-form
        self.F = self.M.diff_form(sg.Integer(2), name="F")
        for i in range(self.n):
            for j in range(self.n):
                self.F[self.ed, i, j] = self.g(self.J(self.ed[i]), self.ed[j])
        self.dF = self.F.derivative()
        print("Kahler:", self.dF == 0)
        # compute torsion form
        self.JdF = self.M.diff_form(sg.Integer(3), name="JdF")
        for i in range(self.n):
            for j in range(i+1, self.n):
                for k in range(j+1, self.n):
                    self.JdF[self.ed, i, j, k] = self.dF(self.J(self.ed[i]), self.J(self.ed[j]), self.J(self.ed[k]))
        print("SKT:", self.is_SKT())


    def is_SKT(self):
        assert self.JdF.degree() == 3
        return self.JdF.derivative() == 0
        

    def nijenhuis(self, X, Y):
        """
        return True if Nijenhuis tensor is 0, and False otherwise
        """
        return X.bracket(Y) + self.J(self.J(X).bracket(Y) + X.bracket(self.J(Y))) - self.J(X).bracket(Y)


    def is_integrable(self):
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
    parser.add_argument("-R", "--curvature", action="store_true", help="Compute the Riemann curvature tensor for the connection specified in '-c'")
    parser.add_argument("-r", "--ricci", action="store_true", help="Compute the Ricci curvature for the connection specified in '-c'")
    parser.add_argument("-s", "--scalar", action="store_true", help="Compute the scalar curvature for the connection specified in '-c'")
    parser.add_argument("-S", "--sectional", action="store_true", help="Compute the sectional curvature for the connection specified in '-c'")
    parser.add_argument("-T", "--torsion", action="store_true", help="Compute the torsion for the connection specified in '-c'")

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
    GEPY = GePy(DIM, FRAME)
    # compute metric for orthonormal basis
    GEPY.compute_metric()
    # if complex manifold, parse complex structure, check integrability, compute fundamental 2-form, 
    if COMPLEX_STRUCTURE:
        GEPY.complex_parsing(COMPLEX_STRUCTURE)
    # compute the Levi-Civita connection
    GEPY.levi_civita_connection()
    # compute Gauduchon connection if specified
    if ARGS.connection == "gauduchon":
        GEPY.gauduchon_connection(ARGS.t)
    # compute curvatures
    NABLA = GEPY.get_connection(ARGS.connection)
    if ARGS.curvature:
        CURVATURE = GEPY.compute_curvature(NABLA)
        print("Riemann curvature tensor:", CURVATURE.display())
    if ARGS.ricci:
        RICCI = GEPY.compute_ricci_curvature(NABLA)
        print("Ricci curvature:", RICCI.display())
    if ARGS.scalar:
        SCALAR = GEPY.compute_scalar_curvature(NABLA)
        print("Scalar curvature:", SCALAR)
    if ARGS.sectional:
        SECTIONAL = GEPY.compute_sectional_curvature(NABLA)
        print("Sectional curvature:")
        for s in SECTIONAL:
            print(s, SECTIONAL[s])
    if ARGS.torsion:
        TORSION = GEPY.compute_torsion(NABLA)
        print("Torsion:", TORSION.display())

    # TODO: implement first and second Ricci curvatures for Gauduchon connection
    # TODO: add more asserts

