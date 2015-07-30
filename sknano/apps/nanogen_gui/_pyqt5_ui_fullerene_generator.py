# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'fullerene_generator.ui'
#
# Created by: PyQt5 UI code generator 5.5
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_FullereneGenerator(object):
    def setupUi(self, FullereneGenerator):
        FullereneGenerator.setObjectName("FullereneGenerator")
        FullereneGenerator.resize(300, 400)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(FullereneGenerator.sizePolicy().hasHeightForWidth())
        FullereneGenerator.setSizePolicy(sizePolicy)
        FullereneGenerator.setMinimumSize(QtCore.QSize(300, 400))
        FullereneGenerator.setMaximumSize(QtCore.QSize(800, 1000))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        FullereneGenerator.setFont(font)
        self.centralwidget = QtWidgets.QWidget(FullereneGenerator)
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
        FullereneGenerator.setCentralWidget(self.centralwidget)

        self.retranslateUi(FullereneGenerator)
        QtCore.QMetaObject.connectSlotsByName(FullereneGenerator)

    def retranslateUi(self, FullereneGenerator):
        _translate = QtCore.QCoreApplication.translate
        FullereneGenerator.setWindowTitle(_translate("FullereneGenerator", "Fullerene Generator"))
        self.label_34.setText(_translate("FullereneGenerator", "Select Fullerene from list..."))

