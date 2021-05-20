# coding: utf8
import Sofa
import SofaRuntime
import numpy as np
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # for running without installing
from SofaOxygenDiffusion import UptakeForceField


class Simulation(object):
    def __init__(self, stl_file='tumor_mesh.msh',
                 dt=10,
                 diffusion_coef=1e6,
                 o2_supply=60,
                 steady_state_threshold=0.005):
        self.stl_file = stl_file
        self.dt = dt
        self.diffusion_coef = diffusion_coef
        self.o2_supply = o2_supply
        self.steady_state_threshold = steady_state_threshold

        self.root_node = Sofa.Core.Node("root")
        self.root_node.addObject("RequiredPlugin", name="SofaOpenglVisual")  # visual stuff
        self.root_node.addObject("RequiredPlugin", name="SofaLoader")  # geometry loaders
        self.root_node.addObject("RequiredPlugin", name="SofaSimpleFem")  # diffusion fem
        self.root_node.addObject("RequiredPlugin", name="SofaBoundaryCondition")  # constraints
        self.root_node.addObject("RequiredPlugin", name="SofaEngine")  # Box Roi
        self.root_node.addObject("RequiredPlugin", name="SofaImplicitOdeSolver")  # implicit solver
        self.root_node.addObject("RequiredPlugin", name="SofaExplicitOdeSolver")  # explicit solver
        self.root_node.addObject("RequiredPlugin", name="SofaMiscForceField")  # meshmatrix
        self.root_node.addObject("RequiredPlugin", name="SofaGeneralEngine")  # TextureInterpolation
        self.root_node.addObject("RequiredPlugin", name="CImgPlugin")  # for loading a bmp image for texture
        self.root_node.addObject("RequiredPlugin", name="SofaBaseLinearSolver")
        self.root_node.addObject("RequiredPlugin", name="SofaGeneralVisual")
        self.root_node.addObject("RequiredPlugin", name="LMConstraint")
        self.root_node.addObject("RequiredPlugin", name="SofaGeneralLoader")

        self.root_node.dt = self.dt
        self.root_node.gravity = [0, 0, 0]
        self.root_node.addObject('VisualStyle', displayFlags="showAll showWireframe")
        self.root_node.addObject('FreeMotionAnimationLoop')
        self.root_node.addObject('GenericConstraintSolver', maxIterations=1000, tolerance=0.0001)

        mesh = self.root_node.addObject('MeshGmshLoader', name='3D_mesh', filename=self.stl_file)
        mesh.position = mesh.position.array()*1e-6
        self.topology = self.root_node.addObject('TetrahedronSetTopologyContainer', name='topology', src='@3D_mesh', tags='mechanics')
        self.topology.init()  # done so we can read how many vertices
        self.root_node.addObject('MechanicalObject', template='Vec3d', name='3D_DOFs', src='@3D_mesh', tags='mechanics')
        self.root_node.addObject('TetrahedronSetGeometryAlgorithms', template='Vec3d', name='GeomAlgo', tags='mechanics')

        oxygen_node = self.root_node.addChild('oxygen')
        oxygen_node.addObject('EulerImplicitSolver', name="Euler", firstOrder=True, tags="oxygen", rayleighStiffness=0.0,
                              rayleighMass=0.0, trapezoidalScheme=False)
        oxygen_node.addObject('CGLinearSolver', name="CG", iterations=1000, tolerance=1.0e-10, threshold=0, tags="oxygen")
        # manually create a mechanical object to initialize values to zero.
        loaded_temps = np.full(self.topology.nbPoints.value, 0)
        self.oxygen = oxygen_node.addObject('MechanicalObject', template="Vec1d", name='oxygen_DOFs',
                                            position=loaded_temps.tolist(), tags="oxygen")

        topology_tetra = self.topology.tetrahedra.array()
        uptake_coef_array = np.full(len(topology_tetra), 0, dtype=np.float32)

        self.uptake_force_field = UptakeForceField(name="UptakeField")
        self.uptake_force_field.set_tetra_uptake_coefficients(uptake_coef_array)
        oxygen_node.addObject(self.uptake_force_field)
        oxygen_node.addObject("TetrahedronDiffusionFEMForceField", name='diffusion', template='Vec1d',
                              constantDiffusionCoefficient=self.diffusion_coef,
                              tagMechanics='mechanics', tags='oxygen')
        oxygen_node.addObject('UncoupledConstraintCorrection')

        oxygen_node.addObject('MeshMatrixMass', template="Vec1d,double", name="Mass", lumping=0, massDensity=1e-2,
                              printLog=False, tags="oxygen")

        visual_node = oxygen_node.addChild('visual_node')
        visual_node.addObject('TextureInterpolation', template="Vec1d", name="EngineInterpolation",
                              input_states="@../oxygen_DOFs.position",
                              input_coordinates="@../../3D_DOFs.position",
                              manual_scale=0, drawPotentiels=True, showIndicesScale=0)
        visual_node.addObject('OglModel', template="Vec3d", name="oglPotential",
                              texcoords="@EngineInterpolation.output_coordinates", texturename="textures/heatColor.bmp",
                              scale3d=[1, 1, 1],
                              material="Default Diffuse 1 1 1 1 0.5 Ambient 1 1 1 1 0.3 Specular 0 0.5 0.5 0.5 1 Emissive\
                              0 0.5 0.5 0.5 1 Shininess 0 45 No texture linked to the material No bump texture linked to\
                              the material ")
        visual_node.addObject('IdentityMapping', input="@../../3D_DOFs", output="@oglPotential")

        self.camera = self.root_node.addObject("InteractiveCamera",
                                               name="camera",
                                               position=[1.21156216e-03,  1.96106548e-04, -4.01305086e-05],
                                               orientation=[-0.08777735,  0.74885405,  0.10808765,  0.64794275],
                                               lookAt=[0, 0, 0],
                                               distance=37,
                                               fieldOfView=45)

        links = self.root_node.addObject('BoxROI', name="box-links", box=[-2.6e-4, -2.6e-4, 2.6e-4, 2.6e-4, 2.6e-4, 2.4e-4],
                                         drawBoxes=True, position="@3D_DOFs.position")

        # set initial temp on left side
        links.init()
        current_02 = self.oxygen.position.array().copy()
        current_02[links.indices.array()] = self.o2_supply
        self.oxygen.position = current_02.tolist()
        oxygen_node.addObject('FixedConstraint', name='flow_in', template="Vec1d", indices="@../box-links.indices")

        for i in range(len(self.oxygen.position.array())):
            oxygen_node.addObject('StopperConstraint', name='oxygen_low'+str(i), template="Vec1d", index=i, min=1e-6,
                                  max=200)

    def initialize(self):
        Sofa.Simulation.init(self.root_node)

    def reset(self):
        Sofa.Simulation.reset(self.root_node)

    def step(self):
        Sofa.Simulation.animate(self.root_node, self.root_node.getDt())
        Sofa.Simulation.updateVisual(self.root_node)

    def solve_steady_state(self):
        """
        Solve system for steady state
        :return: numpy array of tetrahedron oxygen values averaged from each vertex
        """
        prev_sol = self.oxygen.position.array().copy()
        delta = 1e9
        while delta > self.steady_state_threshold:
            Sofa.Simulation.animate(self.root_node, self.root_node.getDt())
            curr_sol = self.oxygen.position.array()
            delta = np.max(np.abs(curr_sol - prev_sol))
            prev_sol = curr_sol.copy()
        Sofa.Simulation.updateVisual(self.root_node)

        oxygen = self.oxygen.position.array()
        tetra = self.topology.tetrahedra.array()
        np_ox = np.zeros([len(tetra), 1])
        for i in range(0, len(tetra)):
            for j in range(0, 4):
                index = tetra[i][j]
                np_ox[i] += 0.25 * oxygen[index]
        return np_ox