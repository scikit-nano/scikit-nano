# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'bulk_structure_generator.ui'
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

class Ui_BulkStructureGenerator(object):
    def setupUi(self, BulkStructureGenerator):
        BulkStructureGenerator.setObjectName(_fromUtf8("BulkStructureGenerator"))
        BulkStructureGenerator.resize(350, 400)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(BulkStructureGenerator.sizePolicy().hasHeightForWidth())
        BulkStructureGenerator.setSizePolicy(sizePolicy)
        BulkStructureGenerator.setMinimumSize(QtCore.QSize(350, 400))
        BulkStructureGenerator.setMaximumSize(QtCore.QSize(800, 1000))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        BulkStructureGenerator.setFont(font)
        self.centralwidget = QtGui.QWidget(BulkStructureGenerator)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout_6 = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.label = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout_3.addWidget(self.label)
        self.bulk_structure_list_widget = QtGui.QListWidget(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.bulk_structure_list_widget.setFont(font)
        self.bulk_structure_list_widget.setObjectName(_fromUtf8("bulk_structure_list_widget"))
        self.verticalLayout_3.addWidget(self.bulk_structure_list_widget)
        self.horizontalLayout_6.addLayout(self.verticalLayout_3)
        self.line = QtGui.QFrame(self.centralwidget)
        self.line.setFrameShape(QtGui.QFrame.VLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.horizontalLayout_6.addWidget(self.line)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_2 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.label_2.setFont(font)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout.addWidget(self.label_2)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.label_3 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(16)
        self.label_3.setFont(font)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.horizontalLayout.addWidget(self.label_3)
        self.nx_spin_box = QtGui.QSpinBox(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.nx_spin_box.setFont(font)
        self.nx_spin_box.setMinimum(1)
        self.nx_spin_box.setMaximum(999)
        self.nx_spin_box.setObjectName(_fromUtf8("nx_spin_box"))
        self.horizontalLayout.addWidget(self.nx_spin_box)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.label_7 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(16)
        self.label_7.setFont(font)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.horizontalLayout_5.addWidget(self.label_7)
        self.ny_spin_box = QtGui.QSpinBox(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.ny_spin_box.setFont(font)
        self.ny_spin_box.setMinimum(1)
        self.ny_spin_box.setMaximum(999)
        self.ny_spin_box.setObjectName(_fromUtf8("ny_spin_box"))
        self.horizontalLayout_5.addWidget(self.ny_spin_box)
        self.verticalLayout.addLayout(self.horizontalLayout_5)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.label_6 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(16)
        self.label_6.setFont(font)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.horizontalLayout_4.addWidget(self.label_6)
        self.nz_spin_box = QtGui.QSpinBox(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        self.nz_spin_box.setFont(font)
        self.nz_spin_box.setMinimum(1)
        self.nz_spin_box.setMaximum(999)
        self.nz_spin_box.setObjectName(_fromUtf8("nz_spin_box"))
        self.horizontalLayout_4.addWidget(self.nz_spin_box)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem)
        self.horizontalLayout_6.addLayout(self.verticalLayout_2)
        BulkStructureGenerator.setCentralWidget(self.centralwidget)

        self.retranslateUi(BulkStructureGenerator)
        QtCore.QMetaObject.connectSlotsByName(BulkStructureGenerator)

    def retranslateUi(self, BulkStructureGenerator):
        BulkStructureGenerator.setWindowTitle(_translate("BulkStructureGenerator", "Bulk Structure Generator", None))
        self.label.setText(_translate("BulkStructureGenerator", "Select Bulk Structure:", None))
        self.label_2.setText(_translate("BulkStructureGenerator", "Unit Cells:", None))
        self.label_3.setText(_translate("BulkStructureGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">n</span><span style=\" vertical-align:sub;\">x </span>=</p></body></html>", None))
        self.label_7.setText(_translate("BulkStructureGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">n</span><span style=\" vertical-align:sub;\">y </span>=</p></body></html>", None))
        self.label_6.setText(_translate("BulkStructureGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">n</span><span style=\" vertical-align:sub;\">z </span>=</p></body></html>", None))

