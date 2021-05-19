import numpy as np
from qtpy.QtWidgets import *
from qtpy.QtGui import *
from qtpy.QtCore import *
from QSofaGLViewTools import QSofaGLView, QSofaViewKeyboardController  # https://github.com/psomers3/QSofaGLViewTools
from tissue_sim import Simulation
from matlab import engine as matlab_engine
import matlab


class MatlabSim(QObject):
    step_finished = Signal(object)  # returns coefficients

    def __init__(self):
        super(MatlabSim, self).__init__()
        self.eng = None

    def start(self):
        self.eng = matlab_engine.start_matlab()
        self.eng.init(nargout=0)

    def run_step(self, oxygen_levels):
        mat_ox = matlab.double(oxygen_levels.tolist())
        self.eng.workspace['oxygen'] = mat_ox
        self.eng.funsim(nargout=0)
        self.step_finished.emit(np.array(self.eng.workspace['tetraUptakeCoefficient']))


class SofaSim(QObject):
    step_finished = Signal(object)  # sends oxygen levels

    def __init__(self):
        super(SofaSim, self).__init__()
        self.sim = Simulation(diffusion_coef=1e6)
        self.forcefield = self.sim.uptake_force_field

    def step_sim(self):
        self.blockSignals(True)
        self.sim.solve_steady_state()
        self.blockSignals(False)

        oxygen = self.sim.oxygen.position.array()
        tetra = self.sim.topology.tetrahedra.array()
        np_ox = np.zeros([len(tetra), 1])
        for i in range(0, len(tetra)):
            for j in range(0, 4):
                index = tetra[i][j]
                np_ox[i] += 0.25 * oxygen[index]
        self.step_finished.emit(np_ox)

    def update_coefficients_and_continue(self, new_coefficients):
        self.forcefield.set_tetra_uptake_coefficients(new_coefficients)
        self.step_sim()


class Simulator(QMainWindow):
    """ Class to coordinate how the simulation runs between the FEM and tumor growth """
    run_sim_step = Signal()

    def __init__(self, autostart=True):
        super(Simulator, self).__init__()

        self.running = False
        self.request_stop = False

        self.sofa_sim = SofaSim()
        self.run_sim_step.connect(self.sofa_sim.step_sim, Qt.QueuedConnection)
        self.matlab_sim = MatlabSim()
        self.matlab_thread = QThread()
        self.matlab_sim.moveToThread(self.matlab_thread)
        self.matlab_thread.start()
        self.matlab_sim.step_finished.connect(self.sofa_sim.update_coefficients_and_continue, Qt.QueuedConnection)
        self.sofa_sim.step_finished.connect(self.matlab_sim.run_step, Qt.QueuedConnection)
        self.matlab_sim.start()

        self.viewer = QSofaGLView(sofa_visuals_node=self.sofa_sim.sim.root_node, camera=self.sofa_sim.sim.camera)
        self.viewer.key_pressed.connect(self.keyPressEvent)
        self.view_control = QSofaViewKeyboardController(translate_rate_limit=0.01)
        self.view_control.set_viewers([self.viewer])
        self.viewer.set_background_color([0, 0, 0, 1])
        self.setCentralWidget(self.viewer)
        # self.view_control.start_auto_update()  # update view continously

        if autostart:
            self.begin()

    def keyPressEvent(self, event: QKeyEvent) -> None:
        key = event.key()
        if key == Qt.Key_Space:
            if self.running:
                self.request_stop = True
            else:
                self.begin()
        elif key == Qt.Key_R:
            self.reset()
        elif key == Qt.Key_P:
            self.save_oxygen_to_file()
        super(Simulator, self).keyPressEvent(event)

    def begin(self):
        self.running = True
        self.run_sim_step.emit()

    def save_oxygen_to_file(self, filename=None):
        if filename is None:
            filename = f"oxygen_at_{self.sofa_sim.sim.root_node.getTime()}"
        np.save(filename, self.sofa_sim.sim.oxygen.position.array())

    def initialize(self):
        self.sofa_sim.sim.initialize()

    def reset(self):
        self.sofa_sim.reset()


if __name__ == '__main__':
    app = QApplication(['sofa_app'])  # viewer engine
    sim = Simulator()
    sim.initialize()
    sim.show()  # show only after sofa sim is initialized
    app.exec()
