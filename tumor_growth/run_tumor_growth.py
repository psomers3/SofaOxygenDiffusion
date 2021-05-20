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
        self.stop_requested = False

    def start(self):
        self.eng = matlab_engine.start_matlab()
        self.eng.init(nargout=0)

    def run_step(self, oxygen_levels):
        if not self.stop_requested:
            mat_ox = matlab.double(oxygen_levels.tolist())
            self.eng.workspace['oxygen'] = mat_ox
            self.eng.funsim(nargout=0)
            self.step_finished.emit(np.array(self.eng.workspace['tetraUptakeCoefficient']))


class SofaSim(QObject):
    step_finished = Signal(object)  # sends oxygen levels

    def __init__(self):
        super(SofaSim, self).__init__()
        self.sim = Simulation(diffusion_coef=1e-4)
        self.forcefield = self.sim.uptake_force_field
        self.stop_requested = False

    def step_sim(self):
        if not self.stop_requested:
            self.blockSignals(True)
            oxygen_levels = self.sim.solve_steady_state()
            self.blockSignals(False)
            self.step_finished.emit(oxygen_levels)

    def update_coefficients_and_continue(self, new_coefficients):
        self.forcefield.set_tetra_uptake_coefficients(new_coefficients)
        self.step_sim()

    def initialize(self):
        self.sim.initialize()


class Simulator(QMainWindow):
    """ Class to coordinate how the simulation runs between the FEM and tumor growth """
    start_sim = Signal()

    def __init__(self, autostart=True):
        super(Simulator, self).__init__()
        self.running = False

        self.sofa_sim = SofaSim()
        self.start_sim.connect(self.sofa_sim.step_sim, Qt.QueuedConnection)

        # create and put tumor growth sim on it's own thread so the GUI isn't locked up when it doesn't need to be.
        self.matlab_sim = MatlabSim()
        self.matlab_thread = QThread()
        self.matlab_sim.moveToThread(self.matlab_thread)
        self.matlab_thread.start()
        self.matlab_sim.step_finished.connect(self.sofa_sim.update_coefficients_and_continue, Qt.QueuedConnection)
        self.sofa_sim.step_finished.connect(self.matlab_sim.run_step, Qt.QueuedConnection)
        self.matlab_sim.start()

        self.viewer = QSofaGLView(sofa_visuals_node=self.sofa_sim.sim.root_node, camera=self.sofa_sim.sim.camera)
        self.view_control = QSofaViewKeyboardController(translate_rate_limit=0.01)
        self.view_control.set_viewers([self.viewer])
        self.viewer.set_background_color([0, 0, 0, 1])
        self.setCentralWidget(self.viewer)
        self.view_control.start_auto_update()  # update view continously

        if autostart:
            self.begin()

    def keyPressEvent(self, event: QKeyEvent) -> None:
        key = event.key()
        if key == Qt.Key_Space:
            if self.running:
                self.stop()
            else:
                self.begin()
        elif key == Qt.Key_R:
            self.reset()
        elif key == Qt.Key_P:
            self.save_oxygen_to_file()
        super(Simulator, self).keyPressEvent(event)

    def begin(self):
        self.running = True
        self.sofa_sim.stop_requested = False
        self.matlab_sim.stop_requested = False
        self.start_sim.emit()

    def stop(self):
        self.running = False
        self.sofa_sim.stop_requested = True
        self.matlab_sim.stop_requested = True

    def save_oxygen_to_file(self, filename=None):
        if filename is None:
            filename = f"oxygen_at_{self.sofa_sim.sim.root_node.getTime()}"
        np.save(filename, self.sofa_sim.sim.oxygen.position.array())

    def initialize(self):
        self.sofa_sim.initialize()

    def reset(self):
        self.sofa_sim.reset()


if __name__ == '__main__':
    app = QApplication(['sofa_app'])
    sim = Simulator(autostart=False)
    sim.initialize()
    sim.show()  # show only after sofa sim is initialized otherwise it will crash
    app.exec()
