# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'nanogen_fullerenes.ui'
#
# Created by: PyQt5 UI code generator 5.4.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_FullereneGeneratorMainWindow(object):
    def setupUi(self, FullereneGeneratorMainWindow):
        FullereneGeneratorMainWindow.setObjectName("FullereneGeneratorMainWindow")
        FullereneGeneratorMainWindow.resize(300, 400)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(FullereneGeneratorMainWindow.sizePolicy().hasHeightForWidth())
        FullereneGeneratorMainWindow.setSizePolicy(sizePolicy)
        FullereneGeneratorMainWindow.setMinimumSize(QtCore.QSize(300, 400))
        FullereneGeneratorMainWindow.setMaximumSize(QtCore.QSize(800, 1000))
        font = QtGui.QFont()
        font.setFamily("Lucida Sans Unicode")
        font.setPointSize(14)
        FullereneGeneratorMainWindow.setFont(font)
        self.centralwidget = QtWidgets.QWidget(FullereneGeneratorMainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.verticalLayout_11 = QtWidgets.QVBoxLayout()
        self.verticalLayout_11.setObjectName("verticalLayout_11")
        self.label_34 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(16)
        self.label_34.setFont(font)
        self.label_34.setObjectName("label_34")
        self.verticalLayout_11.addWidget(self.label_34)
        self.fullerene_list_widget = QtWidgets.QListWidget(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.fullerene_list_widget.setFont(font)
        self.fullerene_list_widget.setObjectName("fullerene_list_widget")
        self.verticalLayout_11.addWidget(self.fullerene_list_widget)
        self.verticalLayout.addLayout(self.verticalLayout_11)
        FullereneGeneratorMainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(FullereneGeneratorMainWindow)
        QtCore.QMetaObject.connectSlotsByName(FullereneGeneratorMainWindow)

    def retranslateUi(self, FullereneGeneratorMainWindow):
        _translate = QtCore.QCoreApplication.translate
        FullereneGeneratorMainWindow.setWindowTitle(_translate("FullereneGeneratorMainWindow", "NanoGen"))
        self.label_34.setText(_translate("FullereneGeneratorMainWindow", "Select Fullerene from list..."))

