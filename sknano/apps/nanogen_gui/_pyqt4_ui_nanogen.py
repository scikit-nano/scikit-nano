# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'nanogen.ui'
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

class Ui_NanoGen(object):
    def setupUi(self, NanoGen):
        NanoGen.setObjectName(_fromUtf8("NanoGen"))
        NanoGen.resize(500, 300)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(NanoGen.sizePolicy().hasHeightForWidth())
        NanoGen.setSizePolicy(sizePolicy)
        NanoGen.setMinimumSize(QtCore.QSize(500, 300))
        NanoGen.setMaximumSize(QtCore.QSize(800, 1000))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        NanoGen.setFont(font)
        self.centralwidget = QtGui.QWidget(NanoGen)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.swnt_generator_push_button = QtGui.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.swnt_generator_push_button.setFont(font)
        self.swnt_generator_push_button.setObjectName(_fromUtf8("swnt_generator_push_button"))
        self.verticalLayout.addWidget(self.swnt_generator_push_button)
        self.mwnt_generator_push_button = QtGui.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.mwnt_generator_push_button.setFont(font)
        self.mwnt_generator_push_button.setObjectName(_fromUtf8("mwnt_generator_push_button"))
        self.verticalLayout.addWidget(self.mwnt_generator_push_button)
        self.graphene_generator_push_button = QtGui.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.graphene_generator_push_button.setFont(font)
        self.graphene_generator_push_button.setObjectName(_fromUtf8("graphene_generator_push_button"))
        self.verticalLayout.addWidget(self.graphene_generator_push_button)
        self.fullerene_generator_push_button = QtGui.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.fullerene_generator_push_button.setFont(font)
        self.fullerene_generator_push_button.setObjectName(_fromUtf8("fullerene_generator_push_button"))
        self.verticalLayout.addWidget(self.fullerene_generator_push_button)
        self.verticalLayout_3.addLayout(self.verticalLayout)
        self.line = QtGui.QFrame(self.centralwidget)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.verticalLayout_3.addWidget(self.line)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.horizontalLayout_25 = QtGui.QHBoxLayout()
        self.horizontalLayout_25.setObjectName(_fromUtf8("horizontalLayout_25"))
        self.save_as_push_button = QtGui.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.save_as_push_button.setFont(font)
        self.save_as_push_button.setObjectName(_fromUtf8("save_as_push_button"))
        self.horizontalLayout_25.addWidget(self.save_as_push_button)
        self.fpath_line_edit = QtGui.QLineEdit(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.fpath_line_edit.setFont(font)
        self.fpath_line_edit.setObjectName(_fromUtf8("fpath_line_edit"))
        self.horizontalLayout_25.addWidget(self.fpath_line_edit)
        self.structure_format_combo_box = QtGui.QComboBox(self.centralwidget)
        self.structure_format_combo_box.setMinimumSize(QtCore.QSize(80, 36))
        self.structure_format_combo_box.setMaximumSize(QtCore.QSize(120, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.structure_format_combo_box.setFont(font)
        self.structure_format_combo_box.setModelColumn(0)
        self.structure_format_combo_box.setObjectName(_fromUtf8("structure_format_combo_box"))
        self.structure_format_combo_box.addItem(_fromUtf8(""))
        self.structure_format_combo_box.addItem(_fromUtf8(""))
        self.horizontalLayout_25.addWidget(self.structure_format_combo_box)
        self.verticalLayout_2.addLayout(self.horizontalLayout_25)
        self.horizontalLayout_15 = QtGui.QHBoxLayout()
        self.horizontalLayout_15.setObjectName(_fromUtf8("horizontalLayout_15"))
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_15.addItem(spacerItem)
        self.generate_push_button = QtGui.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(16)
        self.generate_push_button.setFont(font)
        self.generate_push_button.setObjectName(_fromUtf8("generate_push_button"))
        self.horizontalLayout_15.addWidget(self.generate_push_button)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_15.addItem(spacerItem1)
        self.verticalLayout_2.addLayout(self.horizontalLayout_15)
        self.verticalLayout_3.addLayout(self.verticalLayout_2)
        spacerItem2 = QtGui.QSpacerItem(20, 5, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem2)
        NanoGen.setCentralWidget(self.centralwidget)
        self.nanogen_status_bar = QtGui.QStatusBar(NanoGen)
        self.nanogen_status_bar.setObjectName(_fromUtf8("nanogen_status_bar"))
        NanoGen.setStatusBar(self.nanogen_status_bar)

        self.retranslateUi(NanoGen)
        self.structure_format_combo_box.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(NanoGen)

    def retranslateUi(self, NanoGen):
        NanoGen.setWindowTitle(_translate("NanoGen", "NanoGen", None))
        self.swnt_generator_push_button.setText(_translate("NanoGen", "SWNT Generator", None))
        self.mwnt_generator_push_button.setText(_translate("NanoGen", "MWNT Generator", None))
        self.graphene_generator_push_button.setText(_translate("NanoGen", "Graphene Generator", None))
        self.fullerene_generator_push_button.setText(_translate("NanoGen", "Fullerene Generator", None))
        self.save_as_push_button.setText(_translate("NanoGen", "Save As...", None))
        self.structure_format_combo_box.setItemText(0, _translate("NanoGen", "data", None))
        self.structure_format_combo_box.setItemText(1, _translate("NanoGen", "xyz", None))
        self.generate_push_button.setText(_translate("NanoGen", "Generate", None))

