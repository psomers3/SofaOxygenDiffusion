import numpy as np
from qtpy.QtWidgets import *
from qtpy.QtCore import *
from QSofaGLViewTools import QSofaGLView, QSofaViewKeyboardController
from tissue_sim import Simulation
from matlab import engine as matlab_engine
import matlab


class Simulator(QObject):
    run_matlab_step = Signal()
    run_sim_step = Signal()
    step_completed = Signal()

    def __init__(self, autostart=True):
        super(Simulator, self).__init__()
        self.sim = Simulation(diffusion_coef=1e6)
        self.forcefield = self.sim.uptake_force_field
        self.running = False
        self.request_stop = False
        self.run_matlab_step.connect(self.run_matlab, Qt.QueuedConnection)
        self.run_sim_step.connect(self.step_sim, Qt.QueuedConnection)

        self.eng = matlab_engine.start_matlab()
        self.eng.init(nargout=0)
        if autostart:
            self.begin()

    def run_matlab(self):
        oxygen = self.sim.oxygen.position.array()
        tetra = self.sim.topology.tetrahedra.array()
        np_ox = np.zeros([len(tetra), 1])
        for i in range(0, len(tetra)):
            for j in range(0, 4):
                index = tetra[i][j]
                np_ox[i] += 0.25 * oxygen[index]

        mat_ox = matlab.double(np_ox.tolist())
        self.eng.workspace['oxygen'] = mat_ox
        self.eng.funsim(nargout=0)

        self.matlab_uptake_array = self.eng.workspace['tetraUptakeCoefficient']
        uptake_array = np.array(self.matlab_uptake_array)

        self.forcefield.set_tetra_uptake_coefficients(uptake_array)
        if not self.request_stop:
            self.run_sim_step.emit()
        else:
            self.running = False

    def begin(self):
        self.running = True
        self.run_sim_step.emit()

    def step_sim(self):
        self.now += 1
        self.blockSignals(True)
        self.sim.solve_steady_state()
        self.blockSignals(False)
        self.step_completed.emit()
        self.run_sim_step.emit()

        if not self.request_stop:
            self.run_matlab_step.emit()
        else:
            self.running = False

    def save_oxygen_to_file(self, filename=None):
        if filename is None:
            filename = f"oxygen_at_{self.sim.root_node.getTime()}"
        np.save(filename, self.sim.oxygen.position.array())

    def initialize(self):
        self.sim.initialize()

    def reset(self):
        self.reset()


if __name__ == '__main__':
    app = QApplication(['sofa_app'])  # viewer engine
    sim = Simulator()
    sim.initialize()

    def key_pressed(event):
        key = event.key()
        if key == Qt.Key_Space:
            if sim.running:
                sim.request_stop = True
            else:
                sim.begin()
        elif key == Qt.Key_R:
            sim.reset()
        elif key == Qt.Key_P:
            sim.save_oxygen_to_file()

    viewer = QSofaGLView(sim.sim.root_node, camera=sim.sim.camera)
    view_control = QSofaViewKeyboardController(translate_rate_limit=0.01)
    view_control.set_viewers([viewer])
    # view_control.start_auto_update()  # update view continously, low update rate because sim runs slowly
    sim.step_completed.connect(viewer.update)
    viewer.set_background_color([0, 0, 0, 1])
    viewer.key_pressed.connect(key_pressed)
    viewer.show()

    app.exec()
