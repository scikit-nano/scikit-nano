# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'graphene_generator.ui'
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

class Ui_GrapheneGenerator(object):
    def setupUi(self, GrapheneGenerator):
        GrapheneGenerator.setObjectName(_fromUtf8("GrapheneGenerator"))
        GrapheneGenerator.resize(600, 600)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GrapheneGenerator.sizePolicy().hasHeightForWidth())
        GrapheneGenerator.setSizePolicy(sizePolicy)
        GrapheneGenerator.setMinimumSize(QtCore.QSize(600, 600))
        GrapheneGenerator.setMaximumSize(QtCore.QSize(800, 1000))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        GrapheneGenerator.setFont(font)
        self.centralwidget = QtGui.QWidget(GrapheneGenerator)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayout_7 = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout_7.setObjectName(_fromUtf8("verticalLayout_7"))
        self.verticalLayout_6 = QtGui.QVBoxLayout()
        self.verticalLayout_6.setObjectName(_fromUtf8("verticalLayout_6"))
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.label_33 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_33.setFont(font)
        self.label_33.setObjectName(_fromUtf8("label_33"))
        self.verticalLayout_3.addWidget(self.label_33)
        self.line_10 = QtGui.QFrame(self.centralwidget)
        self.line_10.setFrameShape(QtGui.QFrame.HLine)
        self.line_10.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_10.setObjectName(_fromUtf8("line_10"))
        self.verticalLayout_3.addWidget(self.line_10)
        self.verticalLayout_15 = QtGui.QVBoxLayout()
        self.verticalLayout_15.setObjectName(_fromUtf8("verticalLayout_15"))
        self.conventional_unit_cell_radio_button = QtGui.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.conventional_unit_cell_radio_button.setFont(font)
        self.conventional_unit_cell_radio_button.setChecked(True)
        self.conventional_unit_cell_radio_button.setObjectName(_fromUtf8("conventional_unit_cell_radio_button"))
        self.graphene_generator_button_group = QtGui.QButtonGroup(GrapheneGenerator)
        self.graphene_generator_button_group.setObjectName(_fromUtf8("graphene_generator_button_group"))
        self.graphene_generator_button_group.addButton(self.conventional_unit_cell_radio_button)
        self.verticalLayout_15.addWidget(self.conventional_unit_cell_radio_button)
        self.verticalLayout_13 = QtGui.QVBoxLayout()
        self.verticalLayout_13.setObjectName(_fromUtf8("verticalLayout_13"))
        self.horizontalLayout_19 = QtGui.QHBoxLayout()
        self.horizontalLayout_19.setObjectName(_fromUtf8("horizontalLayout_19"))
        self.horizontalLayout_16 = QtGui.QHBoxLayout()
        self.horizontalLayout_16.setObjectName(_fromUtf8("horizontalLayout_16"))
        self.label_9 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_9.setFont(font)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.horizontalLayout_16.addWidget(self.label_9)
        self.armchair_edge_length_double_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.armchair_edge_length_double_spin_box.setMinimumSize(QtCore.QSize(100, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.armchair_edge_length_double_spin_box.setFont(font)
        self.armchair_edge_length_double_spin_box.setDecimals(4)
        self.armchair_edge_length_double_spin_box.setMinimum(0.1)
        self.armchair_edge_length_double_spin_box.setMaximum(100.0)
        self.armchair_edge_length_double_spin_box.setProperty("value", 1.0)
        self.armchair_edge_length_double_spin_box.setObjectName(_fromUtf8("armchair_edge_length_double_spin_box"))
        self.horizontalLayout_16.addWidget(self.armchair_edge_length_double_spin_box)
        self.horizontalLayout_19.addLayout(self.horizontalLayout_16)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_19.addItem(spacerItem)
        self.verticalLayout_13.addLayout(self.horizontalLayout_19)
        self.horizontalLayout_20 = QtGui.QHBoxLayout()
        self.horizontalLayout_20.setObjectName(_fromUtf8("horizontalLayout_20"))
        self.horizontalLayout_17 = QtGui.QHBoxLayout()
        self.horizontalLayout_17.setObjectName(_fromUtf8("horizontalLayout_17"))
        self.label_10 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_10.setFont(font)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.horizontalLayout_17.addWidget(self.label_10)
        self.zigzag_edge_length_double_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.zigzag_edge_length_double_spin_box.setMinimumSize(QtCore.QSize(100, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.zigzag_edge_length_double_spin_box.setFont(font)
        self.zigzag_edge_length_double_spin_box.setDecimals(4)
        self.zigzag_edge_length_double_spin_box.setMinimum(0.1)
        self.zigzag_edge_length_double_spin_box.setMaximum(100.0)
        self.zigzag_edge_length_double_spin_box.setProperty("value", 1.0)
        self.zigzag_edge_length_double_spin_box.setObjectName(_fromUtf8("zigzag_edge_length_double_spin_box"))
        self.horizontalLayout_17.addWidget(self.zigzag_edge_length_double_spin_box)
        self.horizontalLayout_20.addLayout(self.horizontalLayout_17)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_20.addItem(spacerItem1)
        self.verticalLayout_13.addLayout(self.horizontalLayout_20)
        self.verticalLayout_15.addLayout(self.verticalLayout_13)
        self.verticalLayout_3.addLayout(self.verticalLayout_15)
        self.line_11 = QtGui.QFrame(self.centralwidget)
        self.line_11.setFrameShape(QtGui.QFrame.HLine)
        self.line_11.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_11.setObjectName(_fromUtf8("line_11"))
        self.verticalLayout_3.addWidget(self.line_11)
        self.verticalLayout_17 = QtGui.QVBoxLayout()
        self.verticalLayout_17.setObjectName(_fromUtf8("verticalLayout_17"))
        self.primitive_unit_cell_radio_button = QtGui.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.primitive_unit_cell_radio_button.setFont(font)
        self.primitive_unit_cell_radio_button.setChecked(False)
        self.primitive_unit_cell_radio_button.setObjectName(_fromUtf8("primitive_unit_cell_radio_button"))
        self.graphene_generator_button_group.addButton(self.primitive_unit_cell_radio_button)
        self.verticalLayout_17.addWidget(self.primitive_unit_cell_radio_button)
        self.verticalLayout_19 = QtGui.QVBoxLayout()
        self.verticalLayout_19.setObjectName(_fromUtf8("verticalLayout_19"))
        self.horizontalLayout_23 = QtGui.QHBoxLayout()
        self.horizontalLayout_23.setObjectName(_fromUtf8("horizontalLayout_23"))
        self.horizontalLayout_25 = QtGui.QHBoxLayout()
        self.horizontalLayout_25.setObjectName(_fromUtf8("horizontalLayout_25"))
        self.label_12 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_12.setFont(font)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.horizontalLayout_25.addWidget(self.label_12)
        self.edge_length_double_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.edge_length_double_spin_box.setMinimumSize(QtCore.QSize(100, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.edge_length_double_spin_box.setFont(font)
        self.edge_length_double_spin_box.setReadOnly(True)
        self.edge_length_double_spin_box.setDecimals(4)
        self.edge_length_double_spin_box.setMinimum(0.1)
        self.edge_length_double_spin_box.setMaximum(100.0)
        self.edge_length_double_spin_box.setProperty("value", 1.0)
        self.edge_length_double_spin_box.setObjectName(_fromUtf8("edge_length_double_spin_box"))
        self.horizontalLayout_25.addWidget(self.edge_length_double_spin_box)
        self.horizontalLayout_23.addLayout(self.horizontalLayout_25)
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_23.addItem(spacerItem2)
        self.verticalLayout_19.addLayout(self.horizontalLayout_23)
        self.verticalLayout_17.addLayout(self.verticalLayout_19)
        self.verticalLayout_3.addLayout(self.verticalLayout_17)
        self.line_7 = QtGui.QFrame(self.centralwidget)
        self.line_7.setFrameShape(QtGui.QFrame.HLine)
        self.line_7.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_7.setObjectName(_fromUtf8("line_7"))
        self.verticalLayout_3.addWidget(self.line_7)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.unrolled_swnt_chirality_radio_button = QtGui.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.unrolled_swnt_chirality_radio_button.setFont(font)
        self.unrolled_swnt_chirality_radio_button.setObjectName(_fromUtf8("unrolled_swnt_chirality_radio_button"))
        self.graphene_generator_button_group.addButton(self.unrolled_swnt_chirality_radio_button)
        self.verticalLayout_2.addWidget(self.unrolled_swnt_chirality_radio_button)
        self.horizontalLayout_51 = QtGui.QHBoxLayout()
        self.horizontalLayout_51.setObjectName(_fromUtf8("horizontalLayout_51"))
        self.horizontalLayout_52 = QtGui.QHBoxLayout()
        self.horizontalLayout_52.setObjectName(_fromUtf8("horizontalLayout_52"))
        self.horizontalLayout_53 = QtGui.QHBoxLayout()
        self.horizontalLayout_53.setObjectName(_fromUtf8("horizontalLayout_53"))
        self.label_27 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_27.setFont(font)
        self.label_27.setObjectName(_fromUtf8("label_27"))
        self.horizontalLayout_53.addWidget(self.label_27)
        self.swnt_n_spin_box = QtGui.QSpinBox(self.centralwidget)
        self.swnt_n_spin_box.setMinimumSize(QtCore.QSize(50, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.swnt_n_spin_box.setFont(font)
        self.swnt_n_spin_box.setReadOnly(True)
        self.swnt_n_spin_box.setMaximum(100)
        self.swnt_n_spin_box.setProperty("value", 10)
        self.swnt_n_spin_box.setObjectName(_fromUtf8("swnt_n_spin_box"))
        self.horizontalLayout_53.addWidget(self.swnt_n_spin_box)
        self.horizontalLayout_52.addLayout(self.horizontalLayout_53)
        self.horizontalLayout_54 = QtGui.QHBoxLayout()
        self.horizontalLayout_54.setObjectName(_fromUtf8("horizontalLayout_54"))
        self.label_28 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_28.setFont(font)
        self.label_28.setObjectName(_fromUtf8("label_28"))
        self.horizontalLayout_54.addWidget(self.label_28)
        self.swnt_m_spin_box = QtGui.QSpinBox(self.centralwidget)
        self.swnt_m_spin_box.setMinimumSize(QtCore.QSize(50, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.swnt_m_spin_box.setFont(font)
        self.swnt_m_spin_box.setReadOnly(True)
        self.swnt_m_spin_box.setMaximum(100)
        self.swnt_m_spin_box.setProperty("value", 10)
        self.swnt_m_spin_box.setObjectName(_fromUtf8("swnt_m_spin_box"))
        self.horizontalLayout_54.addWidget(self.swnt_m_spin_box)
        self.horizontalLayout_52.addLayout(self.horizontalLayout_54)
        self.horizontalLayout_51.addLayout(self.horizontalLayout_52)
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_51.addItem(spacerItem3)
        self.verticalLayout_2.addLayout(self.horizontalLayout_51)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.horizontalLayout_58 = QtGui.QHBoxLayout()
        self.horizontalLayout_58.setObjectName(_fromUtf8("horizontalLayout_58"))
        self.label_31 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_31.setFont(font)
        self.label_31.setObjectName(_fromUtf8("label_31"))
        self.horizontalLayout_58.addWidget(self.label_31)
        self.unrolled_swnt_n1_spin_box = QtGui.QSpinBox(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.unrolled_swnt_n1_spin_box.sizePolicy().hasHeightForWidth())
        self.unrolled_swnt_n1_spin_box.setSizePolicy(sizePolicy)
        self.unrolled_swnt_n1_spin_box.setMinimumSize(QtCore.QSize(80, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.unrolled_swnt_n1_spin_box.setFont(font)
        self.unrolled_swnt_n1_spin_box.setReadOnly(True)
        self.unrolled_swnt_n1_spin_box.setMinimum(1)
        self.unrolled_swnt_n1_spin_box.setMaximum(1000)
        self.unrolled_swnt_n1_spin_box.setProperty("value", 1)
        self.unrolled_swnt_n1_spin_box.setObjectName(_fromUtf8("unrolled_swnt_n1_spin_box"))
        self.horizontalLayout_58.addWidget(self.unrolled_swnt_n1_spin_box)
        self.horizontalLayout.addLayout(self.horizontalLayout_58)
        self.horizontalLayout_55 = QtGui.QHBoxLayout()
        self.horizontalLayout_55.setObjectName(_fromUtf8("horizontalLayout_55"))
        self.label_29 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_29.setFont(font)
        self.label_29.setObjectName(_fromUtf8("label_29"))
        self.horizontalLayout_55.addWidget(self.label_29)
        self.unrolled_swnt_n3_spin_box = QtGui.QSpinBox(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.unrolled_swnt_n3_spin_box.sizePolicy().hasHeightForWidth())
        self.unrolled_swnt_n3_spin_box.setSizePolicy(sizePolicy)
        self.unrolled_swnt_n3_spin_box.setMinimumSize(QtCore.QSize(80, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.unrolled_swnt_n3_spin_box.setFont(font)
        self.unrolled_swnt_n3_spin_box.setReadOnly(True)
        self.unrolled_swnt_n3_spin_box.setMinimum(1)
        self.unrolled_swnt_n3_spin_box.setMaximum(1000)
        self.unrolled_swnt_n3_spin_box.setProperty("value", 1)
        self.unrolled_swnt_n3_spin_box.setObjectName(_fromUtf8("unrolled_swnt_n3_spin_box"))
        self.horizontalLayout_55.addWidget(self.unrolled_swnt_n3_spin_box)
        self.horizontalLayout.addLayout(self.horizontalLayout_55)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        spacerItem4 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem4)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.verticalLayout_3.addLayout(self.verticalLayout_2)
        self.horizontalLayout_3.addLayout(self.verticalLayout_3)
        self.line = QtGui.QFrame(self.centralwidget)
        self.line.setFrameShape(QtGui.QFrame.VLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.horizontalLayout_3.addWidget(self.line)
        self.verticalLayout_18 = QtGui.QVBoxLayout()
        self.verticalLayout_18.setObjectName(_fromUtf8("verticalLayout_18"))
        spacerItem5 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_18.addItem(spacerItem5)
        self.verticalLayout_14 = QtGui.QVBoxLayout()
        self.verticalLayout_14.setObjectName(_fromUtf8("verticalLayout_14"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout_21 = QtGui.QHBoxLayout()
        self.horizontalLayout_21.setObjectName(_fromUtf8("horizontalLayout_21"))
        self.horizontalLayout_18 = QtGui.QHBoxLayout()
        self.horizontalLayout_18.setObjectName(_fromUtf8("horizontalLayout_18"))
        self.label_11 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_11.setFont(font)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.horizontalLayout_18.addWidget(self.label_11)
        self.nlayers_spin_box = QtGui.QSpinBox(self.centralwidget)
        self.nlayers_spin_box.setMinimumSize(QtCore.QSize(50, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.nlayers_spin_box.setFont(font)
        self.nlayers_spin_box.setMinimum(1)
        self.nlayers_spin_box.setMaximum(50)
        self.nlayers_spin_box.setObjectName(_fromUtf8("nlayers_spin_box"))
        self.horizontalLayout_18.addWidget(self.nlayers_spin_box)
        self.horizontalLayout_21.addLayout(self.horizontalLayout_18)
        spacerItem6 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_21.addItem(spacerItem6)
        self.verticalLayout.addLayout(self.horizontalLayout_21)
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.label_6 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_6.setFont(font)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.verticalLayout_4.addWidget(self.label_6)
        self.AA_stacking_radio_button = QtGui.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.AA_stacking_radio_button.setFont(font)
        self.AA_stacking_radio_button.setCheckable(True)
        self.AA_stacking_radio_button.setChecked(False)
        self.AA_stacking_radio_button.setObjectName(_fromUtf8("AA_stacking_radio_button"))
        self.stacking_order_button_group = QtGui.QButtonGroup(GrapheneGenerator)
        self.stacking_order_button_group.setObjectName(_fromUtf8("stacking_order_button_group"))
        self.stacking_order_button_group.addButton(self.AA_stacking_radio_button)
        self.verticalLayout_4.addWidget(self.AA_stacking_radio_button)
        self.AB_stacking_radio_button = QtGui.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.AB_stacking_radio_button.setFont(font)
        self.AB_stacking_radio_button.setCheckable(True)
        self.AB_stacking_radio_button.setChecked(True)
        self.AB_stacking_radio_button.setObjectName(_fromUtf8("AB_stacking_radio_button"))
        self.stacking_order_button_group.addButton(self.AB_stacking_radio_button)
        self.verticalLayout_4.addWidget(self.AB_stacking_radio_button)
        self.verticalLayout.addLayout(self.verticalLayout_4)
        self.verticalLayout_14.addLayout(self.verticalLayout)
        self.horizontalLayout_26 = QtGui.QHBoxLayout()
        self.horizontalLayout_26.setObjectName(_fromUtf8("horizontalLayout_26"))
        self.horizontalLayout_22 = QtGui.QHBoxLayout()
        self.horizontalLayout_22.setObjectName(_fromUtf8("horizontalLayout_22"))
        self.label_15 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.label_15.setFont(font)
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.horizontalLayout_22.addWidget(self.label_15)
        self.layer_rotation_increment_double_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.layer_rotation_increment_double_spin_box.setMinimumSize(QtCore.QSize(90, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(14)
        self.layer_rotation_increment_double_spin_box.setFont(font)
        self.layer_rotation_increment_double_spin_box.setReadOnly(True)
        self.layer_rotation_increment_double_spin_box.setMaximum(359.99)
        self.layer_rotation_increment_double_spin_box.setProperty("value", 0.0)
        self.layer_rotation_increment_double_spin_box.setObjectName(_fromUtf8("layer_rotation_increment_double_spin_box"))
        self.horizontalLayout_22.addWidget(self.layer_rotation_increment_double_spin_box)
        self.horizontalLayout_26.addLayout(self.horizontalLayout_22)
        spacerItem7 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_26.addItem(spacerItem7)
        self.verticalLayout_14.addLayout(self.horizontalLayout_26)
        self.verticalLayout_18.addLayout(self.verticalLayout_14)
        spacerItem8 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_18.addItem(spacerItem8)
        self.horizontalLayout_3.addLayout(self.verticalLayout_18)
        self.verticalLayout_6.addLayout(self.horizontalLayout_3)
        self.line_8 = QtGui.QFrame(self.centralwidget)
        self.line_8.setFrameShape(QtGui.QFrame.HLine)
        self.line_8.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_8.setObjectName(_fromUtf8("line_8"))
        self.verticalLayout_6.addWidget(self.line_8)
        self.verticalLayout_5 = QtGui.QVBoxLayout()
        self.verticalLayout_5.setObjectName(_fromUtf8("verticalLayout_5"))
        self.horizontalLayout_31 = QtGui.QHBoxLayout()
        self.horizontalLayout_31.setObjectName(_fromUtf8("horizontalLayout_31"))
        self.horizontalLayout_30 = QtGui.QHBoxLayout()
        self.horizontalLayout_30.setObjectName(_fromUtf8("horizontalLayout_30"))
        self.horizontalLayout_29 = QtGui.QHBoxLayout()
        self.horizontalLayout_29.setObjectName(_fromUtf8("horizontalLayout_29"))
        self.label_17 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(18)
        self.label_17.setFont(font)
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.horizontalLayout_29.addWidget(self.label_17)
        self.element1_combo_box = QtGui.QComboBox(self.centralwidget)
        self.element1_combo_box.setMinimumSize(QtCore.QSize(0, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(18)
        self.element1_combo_box.setFont(font)
        self.element1_combo_box.setObjectName(_fromUtf8("element1_combo_box"))
        self.element1_combo_box.addItem(_fromUtf8(""))
        self.element1_combo_box.addItem(_fromUtf8(""))
        self.element1_combo_box.addItem(_fromUtf8(""))
        self.element1_combo_box.addItem(_fromUtf8(""))
        self.element1_combo_box.addItem(_fromUtf8(""))
        self.horizontalLayout_29.addWidget(self.element1_combo_box)
        self.horizontalLayout_30.addLayout(self.horizontalLayout_29)
        self.horizontalLayout_28 = QtGui.QHBoxLayout()
        self.horizontalLayout_28.setObjectName(_fromUtf8("horizontalLayout_28"))
        self.label_18 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(18)
        self.label_18.setFont(font)
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.horizontalLayout_28.addWidget(self.label_18)
        self.element2_combo_box = QtGui.QComboBox(self.centralwidget)
        self.element2_combo_box.setMinimumSize(QtCore.QSize(0, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(18)
        self.element2_combo_box.setFont(font)
        self.element2_combo_box.setObjectName(_fromUtf8("element2_combo_box"))
        self.element2_combo_box.addItem(_fromUtf8(""))
        self.element2_combo_box.addItem(_fromUtf8(""))
        self.element2_combo_box.addItem(_fromUtf8(""))
        self.element2_combo_box.addItem(_fromUtf8(""))
        self.element2_combo_box.addItem(_fromUtf8(""))
        self.horizontalLayout_28.addWidget(self.element2_combo_box)
        self.horizontalLayout_30.addLayout(self.horizontalLayout_28)
        self.horizontalLayout_31.addLayout(self.horizontalLayout_30)
        spacerItem9 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_31.addItem(spacerItem9)
        self.verticalLayout_5.addLayout(self.horizontalLayout_31)
        self.horizontalLayout_24 = QtGui.QHBoxLayout()
        self.horizontalLayout_24.setObjectName(_fromUtf8("horizontalLayout_24"))
        self.horizontalLayout_14 = QtGui.QHBoxLayout()
        self.horizontalLayout_14.setObjectName(_fromUtf8("horizontalLayout_14"))
        self.elements_bond_label = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(18)
        self.elements_bond_label.setFont(font)
        self.elements_bond_label.setObjectName(_fromUtf8("elements_bond_label"))
        self.horizontalLayout_14.addWidget(self.elements_bond_label)
        self.bond_double_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.bond_double_spin_box.sizePolicy().hasHeightForWidth())
        self.bond_double_spin_box.setSizePolicy(sizePolicy)
        self.bond_double_spin_box.setMinimumSize(QtCore.QSize(150, 36))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(18)
        self.bond_double_spin_box.setFont(font)
        self.bond_double_spin_box.setInputMethodHints(QtCore.Qt.ImhFormattedNumbersOnly)
        self.bond_double_spin_box.setReadOnly(False)
        self.bond_double_spin_box.setDecimals(3)
        self.bond_double_spin_box.setMinimum(0.5)
        self.bond_double_spin_box.setMaximum(10.0)
        self.bond_double_spin_box.setSingleStep(0.01)
        self.bond_double_spin_box.setProperty("value", 1.42)
        self.bond_double_spin_box.setObjectName(_fromUtf8("bond_double_spin_box"))
        self.horizontalLayout_14.addWidget(self.bond_double_spin_box)
        self.horizontalLayout_24.addLayout(self.horizontalLayout_14)
        spacerItem10 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_24.addItem(spacerItem10)
        self.verticalLayout_5.addLayout(self.horizontalLayout_24)
        self.verticalLayout_6.addLayout(self.verticalLayout_5)
        spacerItem11 = QtGui.QSpacerItem(20, 3, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_6.addItem(spacerItem11)
        self.verticalLayout_7.addLayout(self.verticalLayout_6)
        GrapheneGenerator.setCentralWidget(self.centralwidget)

        self.retranslateUi(GrapheneGenerator)
        self.element1_combo_box.setCurrentIndex(3)
        self.element2_combo_box.setCurrentIndex(3)
        QtCore.QMetaObject.connectSlotsByName(GrapheneGenerator)

    def retranslateUi(self, GrapheneGenerator):
        GrapheneGenerator.setWindowTitle(_translate("GrapheneGenerator", "Graphene Generator", None))
        self.label_33.setText(_translate("GrapheneGenerator", "Generate from...", None))
        self.conventional_unit_cell_radio_button.setText(_translate("GrapheneGenerator", "Conventional Unit Cell", None))
        self.label_9.setText(_translate("GrapheneGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">L</span><span style=\" vertical-align:sub;\">armchair</span> =</p></body></html>", None))
        self.armchair_edge_length_double_spin_box.setSuffix(_translate("GrapheneGenerator", " Å", None))
        self.label_10.setText(_translate("GrapheneGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">L</span><span style=\" vertical-align:sub;\">zigzag</span> =</p></body></html>", None))
        self.zigzag_edge_length_double_spin_box.setSuffix(_translate("GrapheneGenerator", " Å", None))
        self.primitive_unit_cell_radio_button.setText(_translate("GrapheneGenerator", "Primitive Unit Cell", None))
        self.label_12.setText(_translate("GrapheneGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">L</span><span style=\" vertical-align:sub;\">edge</span> =</p></body></html>", None))
        self.edge_length_double_spin_box.setSuffix(_translate("GrapheneGenerator", " Å", None))
        self.unrolled_swnt_chirality_radio_button.setText(_translate("GrapheneGenerator", "(n, m) Chiral Indices", None))
        self.label_27.setText(_translate("GrapheneGenerator", "<html><head/><body><p align=\"right\">n = </p></body></html>", None))
        self.label_28.setText(_translate("GrapheneGenerator", "<html><head/><body><p align=\"right\">m = </p></body></html>", None))
        self.label_31.setText(_translate("GrapheneGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">n</span><span style=\" vertical-align:sub;\">1 </span>=</p></body></html>", None))
        self.label_29.setText(_translate("GrapheneGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">n</span><span style=\" vertical-align:sub;\">3 </span>=</p></body></html>", None))
        self.label_11.setText(_translate("GrapheneGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">n</span><span style=\" vertical-align:sub;\">layers</span> =</p></body></html>", None))
        self.label_6.setText(_translate("GrapheneGenerator", "<html><head/><body><p>Stacking Order</p></body></html>", None))
        self.AA_stacking_radio_button.setText(_translate("GrapheneGenerator", "AA", None))
        self.AB_stacking_radio_button.setText(_translate("GrapheneGenerator", "AB", None))
        self.label_15.setText(_translate("GrapheneGenerator", "<html><head/><body><p>deg/layer =</p></body></html>", None))
        self.label_17.setText(_translate("GrapheneGenerator", "Element 1:", None))
        self.element1_combo_box.setItemText(0, _translate("GrapheneGenerator", "H", None))
        self.element1_combo_box.setItemText(1, _translate("GrapheneGenerator", "He", None))
        self.element1_combo_box.setItemText(2, _translate("GrapheneGenerator", "B", None))
        self.element1_combo_box.setItemText(3, _translate("GrapheneGenerator", "C", None))
        self.element1_combo_box.setItemText(4, _translate("GrapheneGenerator", "N", None))
        self.label_18.setText(_translate("GrapheneGenerator", "Element 2:", None))
        self.element2_combo_box.setItemText(0, _translate("GrapheneGenerator", "H", None))
        self.element2_combo_box.setItemText(1, _translate("GrapheneGenerator", "He", None))
        self.element2_combo_box.setItemText(2, _translate("GrapheneGenerator", "B", None))
        self.element2_combo_box.setItemText(3, _translate("GrapheneGenerator", "C", None))
        self.element2_combo_box.setItemText(4, _translate("GrapheneGenerator", "N", None))
        self.elements_bond_label.setText(_translate("GrapheneGenerator", "<html><head/><body><p>C-C bond =</p></body></html>", None))
        self.bond_double_spin_box.setSuffix(_translate("GrapheneGenerator", " Å", None))

