# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mwnt_Ch_list_item_dialog.ui'
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

class Ui_MWNTChListItemDialog(object):
    def setupUi(self, MWNTChListItemDialog):
        MWNTChListItemDialog.setObjectName(_fromUtf8("MWNTChListItemDialog"))
        MWNTChListItemDialog.resize(302, 103)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        MWNTChListItemDialog.setFont(font)
        self.verticalLayout = QtGui.QVBoxLayout(MWNTChListItemDialog)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.label_7 = QtGui.QLabel(MWNTChListItemDialog)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(16)
        self.label_7.setFont(font)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.horizontalLayout_2.addWidget(self.label_7)
        self.n_spin_box = QtGui.QSpinBox(MWNTChListItemDialog)
        self.n_spin_box.setMinimumSize(QtCore.QSize(90, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(16)
        self.n_spin_box.setFont(font)
        self.n_spin_box.setMaximum(100)
        self.n_spin_box.setProperty("value", 10)
        self.n_spin_box.setObjectName(_fromUtf8("n_spin_box"))
        self.horizontalLayout_2.addWidget(self.n_spin_box)
        self.horizontalLayout_3.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.label_14 = QtGui.QLabel(MWNTChListItemDialog)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(16)
        self.label_14.setFont(font)
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.horizontalLayout_4.addWidget(self.label_14)
        self.m_spin_box = QtGui.QSpinBox(MWNTChListItemDialog)
        self.m_spin_box.setMinimumSize(QtCore.QSize(90, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(16)
        self.m_spin_box.setFont(font)
        self.m_spin_box.setMaximum(100)
        self.m_spin_box.setProperty("value", 10)
        self.m_spin_box.setObjectName(_fromUtf8("m_spin_box"))
        self.horizontalLayout_4.addWidget(self.m_spin_box)
        self.horizontalLayout_3.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_5.addLayout(self.horizontalLayout_3)
        self.verticalLayout.addLayout(self.horizontalLayout_5)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.ok_push_button = QtGui.QPushButton(MWNTChListItemDialog)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(12)
        self.ok_push_button.setFont(font)
        self.ok_push_button.setObjectName(_fromUtf8("ok_push_button"))
        self.horizontalLayout.addWidget(self.ok_push_button)
        self.cancel_push_button = QtGui.QPushButton(MWNTChListItemDialog)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(12)
        self.cancel_push_button.setFont(font)
        self.cancel_push_button.setObjectName(_fromUtf8("cancel_push_button"))
        self.horizontalLayout.addWidget(self.cancel_push_button)
        self.verticalLayout.addLayout(self.horizontalLayout)
        spacerItem = QtGui.QSpacerItem(20, 2, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)

        self.retranslateUi(MWNTChListItemDialog)
        QtCore.QMetaObject.connectSlotsByName(MWNTChListItemDialog)

    def retranslateUi(self, MWNTChListItemDialog):
        MWNTChListItemDialog.setWindowTitle(_translate("MWNTChListItemDialog", "(n, m) Dialog", None))
        self.label_7.setText(_translate("MWNTChListItemDialog", "<html><head/><body><p align=\"right\">n = </p></body></html>", None))
        self.label_14.setText(_translate("MWNTChListItemDialog", "<html><head/><body><p align=\"right\">m = </p></body></html>", None))
        self.ok_push_button.setText(_translate("MWNTChListItemDialog", "OK", None))
        self.cancel_push_button.setText(_translate("MWNTChListItemDialog", "Cancel", None))

