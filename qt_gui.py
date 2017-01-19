from PyQt4 import QtGui, QtCore
import SelLabel
import sys


class MainWindow(QtGui.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.setGeometry(50, 50, 800, 500)
        self.setWindowTitle("Combinatorial Labeling")
        self.add_stuff()
        self.show()

    def add_stuff(self):
        button = QtGui.QPushButton("Some button", self)
        button.clicked.connect(self.button_click)
        button.move(50, 50)
        button.adjustSize()

    def button_click(self):
        print("click")

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec())
