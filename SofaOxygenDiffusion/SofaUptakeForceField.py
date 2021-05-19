# coding: utf8
import Sofa
import numpy as np
from numba import jit


@jit(nopython=True)
def _addForce(f, x, x_prior, o2_coef, edge_uptake_coefficients, edges):
    for i in range(0, len(edges)):
        v0 = edges[i][0]
        v1 = edges[i][1]
        c_prior_average = (x_prior[v1] + x_prior[v0]) / 2
        x_average = (x[v1] + x[v0]) / 2  # average the values on each vertex
        dp2 = (edge_uptake_coefficients[i] / (o2_coef + c_prior_average)) * x_average
        f[v1] += dp2
        f[v0] += dp2
    return f


@jit(nopython=True)
def _addDForce(df, dx, x_prior, k_factor, o2_coef, edge_uptake_coefficients, edges):
    for i in range(len(edges)):
        v0 = edges[i][0]
        v1 = edges[i][1]
        dx_average = (dx[v1] + dx[v0]) / 2  # average the values on each vertex
        c_average = (x_prior[v1] + x_prior[v0]) / 2
        dp2 = (edge_uptake_coefficients[i] / (o2_coef + c_average)) * dx_average * k_factor
        df[v1] += dp2
        df[v0] += dp2
    return df


class UptakeForceField(Sofa.Core.ForceFieldVec1d):
    """Implementation of Oxygen uptake in python"""

    def __init__(self, *args, **kwargs):
        Sofa.Core.ForceFieldVec1d.__init__(self, *args, **kwargs)
        self.d_O2SaturationConstant = 2.5  # in mmHG
        self.sim_time = 0

    def init(self):
        self.mechanicalState = self.getContext().getMechanicalState()
        self.c_prior = self.mechanicalState.position.array().copy()
        self.topology = self.getContext().getMeshTopology()
        self.edges = self.topology.edges.array()
        self.nbEdges = len(self.edges)
        self.edgeUptakeCoefficient = np.zeros(self.nbEdges, dtype=np.float32)
        self.tetra = self.topology.tetrahedra.array()
        self.sim_time = 0
        if self.d_tetra_uptake_coefs is None:
            self.set_tetra_uptake_coefficients(np.zeros(len(self.tetra), dtype=float))
        self.computeEdgeDiffusionCoefficient()

    def set_tetra_uptake_coefficients(self, coefs):
        self.d_tetra_uptake_coefs = coefs
        try:
            self.computeEdgeDiffusionCoefficient()
        except KeyError:
            # not initialized yet
            pass

    def computeEdgeDiffusionCoefficient(self):
        nbTetra = len(self.tetra)
        position = self.topology.position.array()
        point = np.zeros((4, 3))
        shapeVector = np.zeros((4, 3))
        local_edges_in_tetrahedron = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])

        for i in range(0, nbTetra):
            tetrahedron = self.tetra[i]
            for j in range(0, 4):
                point[j] = position[tetrahedron[j]]

            # compute 6 times the rest volume
            volume = np.dot(np.cross(point[1] - point[0], point[2] - point[0]), point[0] - point[3])

            # store shape vectors
            for j in range(0, 4):
                if (j % 2) == 0:
                    shapeVector[j] = -np.cross(point[(j + 2) % 4] - point[(j + 1) % 4],
                                               point[(j + 3) % 4] - point[(j + 1) % 4]) / volume
                else:
                    shapeVector[j] = np.cross(point[(j + 2) % 4] - point[(j + 1) % 4],
                                              point[(j + 3) % 4] - point[(j + 1) % 4]) / volume

            diff2 = self.d_tetra_uptake_coefs[i] * np.absolute(volume) / 6

            edges_in_tetrahedron = self.topology.getEdgesInTetrahedron(i)
            for j in range(0, 6):
                k = local_edges_in_tetrahedron[j][0]
                l = local_edges_in_tetrahedron[j][1]
                dot_product = np.dot(shapeVector[k], shapeVector[l])
                val2 = dot_product * diff2
                self.edgeUptakeCoefficient[edges_in_tetrahedron[j]] += val2

    def addForce(self, m, force, c, vel):
        if self.sim_time != self.getContext().getTime():
            self.c_prior = self.mechanicalState.position.array().copy()
            self.sim_time = self.getContext().getTime()
        with force.writeableArray() as f:
            return _addForce(f, c.value, self.c_prior, self.d_O2SaturationConstant, self.edgeUptakeCoefficient,
                             self.edges)

    def addDForce(self, params, dforce, dc):
        if self.sim_time != self.getContext().getTime():
            self.c_prior = self.mechanicalState.position.array().copy()
            self.sim_time = self.getContext().getTime()
        with dforce.writeableArray() as df:
            return _addDForce(df, dc.value, self.c_prior, params['kFactor'], self.d_O2SaturationConstant,
                              self.edgeUptakeCoefficient, self.edges)
