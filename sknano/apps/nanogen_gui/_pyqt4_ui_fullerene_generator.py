# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'fullerene_generator.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_FullereneGenerator(object):
    def setupUi(self, FullereneGenerator):
        FullereneGenerator.setObjectName(_fromUtf8("FullereneGenerator"))
        FullereneGenerator.resize(300, 400)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(FullereneGenerator.sizePolicy().hasHeightForWidth())
        FullereneGenerator.setSizePolicy(sizePolicy)
        FullereneGenerator.setMinimumSize(QtCore.QSize(300, 400))
        FullereneGenerator.setMaximumSize(QtCore.QSize(800, 1000))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        FullereneGenerator.setFont(font)
        self.centralwidget = QtGui.QWidget(FullereneGenerator)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.verticalLayout_11 = QtGui.QVBoxLayout()
        self.verticalLayout_11.setObjectName(_fromUtf8("verticalLayout_11"))
        self.label_34 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(16)
        self.label_34.setFont(font)
        self.label_34.setObjectName(_fromUtf8("label_34"))
        self.verticalLayout_11.addWidget(self.label_34)
        self.fullerene_list_widget = QtGui.QListWidget(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(12)
        self.fullerene_list_widget.setFont(font)
        self.fullerene_list_widget.setObjectName(_fromUtf8("fullerene_list_widget"))
        self.verticalLayout_11.addWidget(self.fullerene_list_widget)
        self.verticalLayout.addLayout(self.verticalLayout_11)
        FullereneGenerator.setCentralWidget(self.centralwidget)

        self.retranslateUi(FullereneGenerator)
        QtCore.QMetaObject.connectSlotsByName(FullereneGenerator)

    def retranslateUi(self, FullereneGenerator):
        FullereneGenerator.setWindowTitle(_translate("FullereneGenerator", "Fullerene Generator", None))
        self.label_34.setText(_translate("FullereneGenerator", "Select Fullerene from list...", None))

