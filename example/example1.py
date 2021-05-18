# coding: utf8
import Sofa
import SofaRuntime
import numpy as np
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from SofaOxygenDiffusion.SofaUptakeForceField import UptakeForceField


if __name__ == '__main__':
    from qtpy.QtWidgets import *
    from qtpy.QtCore import *
    from QSofaGLViewTools import QSofaGLView, QSofaViewKeyboardController


def createScene(node):
    node.addObject("RequiredPlugin", name="SofaOpenglVisual")  # visual stuff
    node.addObject("RequiredPlugin", name="SofaLoader")  # geometry loaders
    node.addObject("RequiredPlugin", name="SofaSimpleFem")  # diffusion fem
    node.addObject("RequiredPlugin", name="SofaBoundaryCondition")  # constraints
    node.addObject("RequiredPlugin", name="SofaEngine")  # Box Roi
    node.addObject("RequiredPlugin", name="SofaImplicitOdeSolver")  # implicit solver
    node.addObject("RequiredPlugin", name="SofaExplicitOdeSolver")  # explicit solver
    node.addObject("RequiredPlugin", name="SofaMiscForceField")  # meshmatrix
    node.addObject("RequiredPlugin", name="SofaGeneralEngine")  # TextureInterpolation
    node.addObject("RequiredPlugin", name="CImgPlugin")  # for loading a bmp image for texture
    node.addObject("RequiredPlugin", name="SofaBaseLinearSolver")
    node.addObject("RequiredPlugin", name="SofaGeneralVisual")
    node.addObject("RequiredPlugin", name="LMConstraint")
    node.addObject("RequiredPlugin", name="SofaGeneralLoader")

    node.dt = 10  # such a low diffusion coefficient makes it take forever to diffuse
    node.gravity = [0, 0, 0]
    node.addObject('VisualStyle', displayFlags="showAll showWireframe")
    node.addObject('FreeMotionAnimationLoop')
    node.addObject('GenericConstraintSolver', maxIterations=1000, tolerance=0.0001)

    mesh = node.addObject('MeshGmshLoader', name='3D_mesh', filename='TetraMesh.msh')
    mesh.position = mesh.position.array()*1e-6
    topology = node.addObject('TetrahedronSetTopologyContainer', name='topology', src='@3D_mesh', tags='mechanics')
    topology.init()  # done so we can read how many vertices
    node.addObject('MechanicalObject', template='Vec3d', name='3D_DOFs', src='@3D_mesh', tags='mechanics')
    node.addObject('TetrahedronSetGeometryAlgorithms', template='Vec3d', name='GeomAlgo', tags='mechanics')

    oxygen_node = node.addChild('oxygen')
    oxygen_node.addObject('EulerImplicitSolver', name="Euler", firstOrder=True, tags="oxygen", rayleighStiffness=0.0,
                          rayleighMass=0.0, trapezoidalScheme=False)
    oxygen_node.addObject('CGLinearSolver', name="CG", iterations=1000, tolerance=1.0e-10, threshold=0, tags="oxygen")
    # manually create a mechanical object to initialize values to zero.
    loaded_02 = np.full(topology.nbPoints.value, 0)
    oxygen = oxygen_node.addObject('MechanicalObject', template="Vec1d", name='oxygen_DOFs',
                                   position=loaded_02.tolist(), tags="oxygen")
    tumor_roi = node.addObject('BoxROI', name="tumor", box=[-1e-4, -1e-4, -1e-4, 1e-4, 1e-4, 1e-4], drawBoxes=True,
                               position="@3D_DOFs.position", tetrahedra="@topology.tetrahedra")
    tumor_roi.init()
    topology_tetra = topology.tetrahedra.array()
    tumor_indices = tumor_roi.tetrahedronIndices.array()

    uptake_coef_array = np.full(len(topology_tetra), 0, dtype=np.float32)
    uptake_coef_array[tumor_indices] = np.full(len(tumor_indices), 0.5, dtype=np.float32)

    force_field = UptakeForceField(name="UptakeField")
    force_field.set_tetra_uptake_coefficients(uptake_coef_array)
    oxygen_node.addObject(force_field)
    oxygen_node.addObject("TetrahedronDiffusionFEMForceField", name='diffusion', template='Vec1d',
                          tagMechanics='mechanics', tags='oxygen')
    oxygen_node.addObject('UncoupledConstraintCorrection')

    oxygen_node.addObject('MeshMatrixMass', template="Vec1d,double", name="Mass", lumping=0, massDensity=1e-12,
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

    node.addObject("InteractiveCamera",
                   name="camera",
                   position=[1.21156216e-03,  1.96106548e-04, -4.01305086e-05],
                   orientation=[-0.08777735,  0.74885405,  0.10808765,  0.64794275],
                   lookAt=[0, 0, 0],
                   distance=37,
                   fieldOfView=45)

    links = node.addObject('BoxROI', name="box-links", box=[-2.6e-4, -2.6e-4, 2.6e-4, 2.6e-4, 2.6e-4, 2.4e-4],
                           drawBoxes=True, position="@3D_DOFs.position")

    # set initial temp on left side
    links.init()
    current_temp = oxygen.position.array().copy()
    current_temp[links.indices.array()] = 2.5
    oxygen.position = current_temp.tolist()
    oxygen_node.addObject('FixedConstraint', name='flow_in', template="Vec1d", indices="@../box-links.indices")

    for i in range(len(oxygen.position.array())):
        oxygen_node.addObject('StopperConstraint', name='oxygen_low'+str(i), template="Vec1d", index=i, min=1e-6,
                              max=200)
    return node


if __name__ == '__main__':
    class SimAnimator(QObject):
        step_finished = Signal()

        def __init__(self, root_node: Sofa.Core.Node):
            super(SimAnimator, self).__init__()
            self.node = root_node
            self.step_finished.connect(self.step_sim, Qt.QueuedConnection)
            self.running = False

        def step_sim(self):
            if not self.running:
                return
            self.blockSignals(True)
            Sofa.Simulation.animate(self.node, self.node.getDt())
            Sofa.Simulation.updateVisual(self.node)
            self.blockSignals(False)
            self.step_finished.emit()

    app = QApplication(['sofa_app'])  # viewer engine

    root_node = Sofa.Core.Node("root")
    createScene(root_node)  # create sofa scene
    Sofa.Simulation.init(root_node)
    viewer = QSofaGLView(root_node, camera=root_node.camera)

    view_control = QSofaViewKeyboardController(translate_rate_limit=0.01)
    view_control.set_viewers([viewer])
    view_control.start_auto_update()  # update view continously
    viewer.show()

    sim = SimAnimator(root_node)

    def key_pressed(event):
        key = event.key()
        if key == Qt.Key_Space:
            if sim.running:
                sim.running = False
            else:
                sim.running = True
                sim.step_finished.emit()
        elif key == Qt.Key_R:
            Sofa.Simulation.reset(root_node)

    viewer.set_background_color([0, 0, 0, 1])
    viewer.key_pressed.connect(key_pressed)
    app.exec()
