# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mwnt_generator.ui'
#
# Created: Thu May  5 11:45:35 2016
#      by: PyQt5 UI code generator 5.2.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MWNTGenerator(object):
    def setupUi(self, MWNTGenerator):
        MWNTGenerator.setObjectName("MWNTGenerator")
        MWNTGenerator.resize(650, 650)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MWNTGenerator.sizePolicy().hasHeightForWidth())
        MWNTGenerator.setSizePolicy(sizePolicy)
        MWNTGenerator.setMinimumSize(QtCore.QSize(650, 650))
        MWNTGenerator.setMaximumSize(QtCore.QSize(800, 1000))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        MWNTGenerator.setFont(font)
        self.centralwidget = QtWidgets.QWidget(MWNTGenerator)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.label_19 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label_19.setFont(font)
        self.label_19.setObjectName("label_19")
        self.verticalLayout_4.addWidget(self.label_19)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.mwnt_Ch_list_radio_button = QtWidgets.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.mwnt_Ch_list_radio_button.setFont(font)
        self.mwnt_Ch_list_radio_button.setChecked(True)
        self.mwnt_Ch_list_radio_button.setObjectName("mwnt_Ch_list_radio_button")
        self.mwnt_generator_button_group = QtWidgets.QButtonGroup(MWNTGenerator)
        self.mwnt_generator_button_group.setObjectName("mwnt_generator_button_group")
        self.mwnt_generator_button_group.addButton(self.mwnt_Ch_list_radio_button)
        self.verticalLayout.addWidget(self.mwnt_Ch_list_radio_button)
        self.mwnt_Ch_list_widget = QtWidgets.QListWidget(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.mwnt_Ch_list_widget.setFont(font)
        self.mwnt_Ch_list_widget.setObjectName("mwnt_Ch_list_widget")
        self.verticalLayout.addWidget(self.mwnt_Ch_list_widget)
        self.add_Ch_push_button = QtWidgets.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.add_Ch_push_button.setFont(font)
        self.add_Ch_push_button.setObjectName("add_Ch_push_button")
        self.verticalLayout.addWidget(self.add_Ch_push_button)
        self.edit_selected_Ch_push_button = QtWidgets.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.edit_selected_Ch_push_button.setFont(font)
        self.edit_selected_Ch_push_button.setObjectName("edit_selected_Ch_push_button")
        self.verticalLayout.addWidget(self.edit_selected_Ch_push_button)
        self.remove_selected_Ch_push_button = QtWidgets.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.remove_selected_Ch_push_button.setFont(font)
        self.remove_selected_Ch_push_button.setObjectName("remove_selected_Ch_push_button")
        self.verticalLayout.addWidget(self.remove_selected_Ch_push_button)
        self.clear_Ch_list_push_button = QtWidgets.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.clear_Ch_list_push_button.setFont(font)
        self.clear_Ch_list_push_button.setObjectName("clear_Ch_list_push_button")
        self.verticalLayout.addWidget(self.clear_Ch_list_push_button)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.line_4 = QtWidgets.QFrame(self.centralwidget)
        self.line_4.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.horizontalLayout.addWidget(self.line_4)
        self.verticalLayout_9 = QtWidgets.QVBoxLayout()
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.mwnt_wall_parameters_radio_button = QtWidgets.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.mwnt_wall_parameters_radio_button.setFont(font)
        self.mwnt_wall_parameters_radio_button.setObjectName("mwnt_wall_parameters_radio_button")
        self.mwnt_generator_button_group.addButton(self.mwnt_wall_parameters_radio_button)
        self.verticalLayout_9.addWidget(self.mwnt_wall_parameters_radio_button)
        self.horizontalLayout_49 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_49.setObjectName("horizontalLayout_49")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout()
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.horizontalLayout_48 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_48.setObjectName("horizontalLayout_48")
        self.horizontalLayout_47 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_47.setObjectName("horizontalLayout_47")
        self.label_24 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.label_24.setFont(font)
        self.label_24.setObjectName("label_24")
        self.horizontalLayout_47.addWidget(self.label_24)
        self.Nwalls_spin_box = QtWidgets.QSpinBox(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.Nwalls_spin_box.setFont(font)
        self.Nwalls_spin_box.setMinimum(1)
        self.Nwalls_spin_box.setProperty("value", 3)
        self.Nwalls_spin_box.setObjectName("Nwalls_spin_box")
        self.horizontalLayout_47.addWidget(self.Nwalls_spin_box)
        self.horizontalLayout_48.addLayout(self.horizontalLayout_47)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_48.addItem(spacerItem)
        self.verticalLayout_8.addLayout(self.horizontalLayout_48)
        self.horizontalLayout_43 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_43.setObjectName("horizontalLayout_43")
        self.horizontalLayout_44 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_44.setObjectName("horizontalLayout_44")
        self.label_25 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.label_25.setFont(font)
        self.label_25.setObjectName("label_25")
        self.horizontalLayout_44.addWidget(self.label_25)
        self.min_wall_diameter_double_spin_box = QtWidgets.QDoubleSpinBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.min_wall_diameter_double_spin_box.sizePolicy().hasHeightForWidth())
        self.min_wall_diameter_double_spin_box.setSizePolicy(sizePolicy)
        self.min_wall_diameter_double_spin_box.setMinimumSize(QtCore.QSize(100, 36))
        self.min_wall_diameter_double_spin_box.setMaximumSize(QtCore.QSize(100, 36))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.min_wall_diameter_double_spin_box.setFont(font)
        self.min_wall_diameter_double_spin_box.setDecimals(1)
        self.min_wall_diameter_double_spin_box.setMinimum(0.0)
        self.min_wall_diameter_double_spin_box.setMaximum(1000.0)
        self.min_wall_diameter_double_spin_box.setProperty("value", 0.0)
        self.min_wall_diameter_double_spin_box.setObjectName("min_wall_diameter_double_spin_box")
        self.horizontalLayout_44.addWidget(self.min_wall_diameter_double_spin_box)
        self.horizontalLayout_43.addLayout(self.horizontalLayout_44)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_43.addItem(spacerItem1)
        self.verticalLayout_8.addLayout(self.horizontalLayout_43)
        self.horizontalLayout_40 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_40.setObjectName("horizontalLayout_40")
        self.horizontalLayout_41 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_41.setObjectName("horizontalLayout_41")
        self.label_23 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.label_23.setFont(font)
        self.label_23.setObjectName("label_23")
        self.horizontalLayout_41.addWidget(self.label_23)
        self.max_wall_diameter_double_spin_box = QtWidgets.QDoubleSpinBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.max_wall_diameter_double_spin_box.sizePolicy().hasHeightForWidth())
        self.max_wall_diameter_double_spin_box.setSizePolicy(sizePolicy)
        self.max_wall_diameter_double_spin_box.setMinimumSize(QtCore.QSize(100, 36))
        self.max_wall_diameter_double_spin_box.setMaximumSize(QtCore.QSize(100, 36))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.max_wall_diameter_double_spin_box.setFont(font)
        self.max_wall_diameter_double_spin_box.setDecimals(1)
        self.max_wall_diameter_double_spin_box.setMinimum(0.1)
        self.max_wall_diameter_double_spin_box.setMaximum(9999.9)
        self.max_wall_diameter_double_spin_box.setProperty("value", 100.0)
        self.max_wall_diameter_double_spin_box.setObjectName("max_wall_diameter_double_spin_box")
        self.horizontalLayout_41.addWidget(self.max_wall_diameter_double_spin_box)
        self.horizontalLayout_40.addLayout(self.horizontalLayout_41)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_40.addItem(spacerItem2)
        self.verticalLayout_8.addLayout(self.horizontalLayout_40)
        self.horizontalLayout_45 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_45.setObjectName("horizontalLayout_45")
        self.horizontalLayout_46 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_46.setObjectName("horizontalLayout_46")
        self.label_26 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.label_26.setFont(font)
        self.label_26.setObjectName("label_26")
        self.horizontalLayout_46.addWidget(self.label_26)
        self.wall_spacing_double_spin_box = QtWidgets.QDoubleSpinBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.wall_spacing_double_spin_box.sizePolicy().hasHeightForWidth())
        self.wall_spacing_double_spin_box.setSizePolicy(sizePolicy)
        self.wall_spacing_double_spin_box.setMinimumSize(QtCore.QSize(100, 36))
        self.wall_spacing_double_spin_box.setMaximumSize(QtCore.QSize(100, 36))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.wall_spacing_double_spin_box.setFont(font)
        self.wall_spacing_double_spin_box.setDecimals(2)
        self.wall_spacing_double_spin_box.setMinimum(0.1)
        self.wall_spacing_double_spin_box.setMaximum(100.0)
        self.wall_spacing_double_spin_box.setProperty("value", 3.35)
        self.wall_spacing_double_spin_box.setObjectName("wall_spacing_double_spin_box")
        self.horizontalLayout_46.addWidget(self.wall_spacing_double_spin_box)
        self.horizontalLayout_45.addLayout(self.horizontalLayout_46)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_45.addItem(spacerItem3)
        self.verticalLayout_8.addLayout(self.horizontalLayout_45)
        spacerItem4 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_8.addItem(spacerItem4)
        self.horizontalLayout_49.addLayout(self.verticalLayout_8)
        self.verticalLayout_7 = QtWidgets.QVBoxLayout()
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_7.addWidget(self.label_5)
        self.verticalLayout_6 = QtWidgets.QVBoxLayout()
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.random_wall_chiralities_radio_button = QtWidgets.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.random_wall_chiralities_radio_button.setFont(font)
        self.random_wall_chiralities_radio_button.setChecked(True)
        self.random_wall_chiralities_radio_button.setObjectName("random_wall_chiralities_radio_button")
        self.mwnts_wall_chiral_type_button_group = QtWidgets.QButtonGroup(MWNTGenerator)
        self.mwnts_wall_chiral_type_button_group.setObjectName("mwnts_wall_chiral_type_button_group")
        self.mwnts_wall_chiral_type_button_group.addButton(self.random_wall_chiralities_radio_button)
        self.verticalLayout_6.addWidget(self.random_wall_chiralities_radio_button)
        self.armchair_wall_chiralities_radio_button = QtWidgets.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.armchair_wall_chiralities_radio_button.setFont(font)
        self.armchair_wall_chiralities_radio_button.setObjectName("armchair_wall_chiralities_radio_button")
        self.mwnts_wall_chiral_type_button_group.addButton(self.armchair_wall_chiralities_radio_button)
        self.verticalLayout_6.addWidget(self.armchair_wall_chiralities_radio_button)
        self.zigzag_wall_chiralities_radio_button = QtWidgets.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.zigzag_wall_chiralities_radio_button.setFont(font)
        self.zigzag_wall_chiralities_radio_button.setObjectName("zigzag_wall_chiralities_radio_button")
        self.mwnts_wall_chiral_type_button_group.addButton(self.zigzag_wall_chiralities_radio_button)
        self.verticalLayout_6.addWidget(self.zigzag_wall_chiralities_radio_button)
        self.achiral_wall_chiralities_radio_button = QtWidgets.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.achiral_wall_chiralities_radio_button.setFont(font)
        self.achiral_wall_chiralities_radio_button.setObjectName("achiral_wall_chiralities_radio_button")
        self.mwnts_wall_chiral_type_button_group.addButton(self.achiral_wall_chiralities_radio_button)
        self.verticalLayout_6.addWidget(self.achiral_wall_chiralities_radio_button)
        self.chiral_wall_chiralities_radio_button = QtWidgets.QRadioButton(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.chiral_wall_chiralities_radio_button.setFont(font)
        self.chiral_wall_chiralities_radio_button.setObjectName("chiral_wall_chiralities_radio_button")
        self.mwnts_wall_chiral_type_button_group.addButton(self.chiral_wall_chiralities_radio_button)
        self.verticalLayout_6.addWidget(self.chiral_wall_chiralities_radio_button)
        self.verticalLayout_7.addLayout(self.verticalLayout_6)
        spacerItem5 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_7.addItem(spacerItem5)
        self.horizontalLayout_49.addLayout(self.verticalLayout_7)
        self.verticalLayout_9.addLayout(self.horizontalLayout_49)
        self.horizontalLayout.addLayout(self.verticalLayout_9)
        self.verticalLayout_4.addLayout(self.horizontalLayout)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.horizontalLayout_34 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_34.setObjectName("horizontalLayout_34")
        self.horizontalLayout_36 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_36.setObjectName("horizontalLayout_36")
        self.label_20 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label_20.setFont(font)
        self.label_20.setObjectName("label_20")
        self.horizontalLayout_36.addWidget(self.label_20)
        self.mwnt_L_double_spin_box = QtWidgets.QDoubleSpinBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mwnt_L_double_spin_box.sizePolicy().hasHeightForWidth())
        self.mwnt_L_double_spin_box.setSizePolicy(sizePolicy)
        self.mwnt_L_double_spin_box.setMinimumSize(QtCore.QSize(150, 36))
        self.mwnt_L_double_spin_box.setMaximumSize(QtCore.QSize(200, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.mwnt_L_double_spin_box.setFont(font)
        self.mwnt_L_double_spin_box.setDecimals(4)
        self.mwnt_L_double_spin_box.setMinimum(0.1)
        self.mwnt_L_double_spin_box.setMaximum(100.0)
        self.mwnt_L_double_spin_box.setProperty("value", 1.0)
        self.mwnt_L_double_spin_box.setObjectName("mwnt_L_double_spin_box")
        self.horizontalLayout_36.addWidget(self.mwnt_L_double_spin_box)
        self.horizontalLayout_34.addLayout(self.horizontalLayout_36)
        spacerItem6 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_34.addItem(spacerItem6)
        self.verticalLayout_3.addLayout(self.horizontalLayout_34)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.bundle_generator_check_box = QtWidgets.QCheckBox(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.bundle_generator_check_box.setFont(font)
        self.bundle_generator_check_box.setCheckable(True)
        self.bundle_generator_check_box.setObjectName("bundle_generator_check_box")
        self.verticalLayout_2.addWidget(self.bundle_generator_check_box)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.horizontalLayout_32 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_32.setObjectName("horizontalLayout_32")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.horizontalLayout_32.addWidget(self.label_8)
        self.bundle_n1_spin_box = QtWidgets.QSpinBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.bundle_n1_spin_box.sizePolicy().hasHeightForWidth())
        self.bundle_n1_spin_box.setSizePolicy(sizePolicy)
        self.bundle_n1_spin_box.setMinimumSize(QtCore.QSize(60, 36))
        self.bundle_n1_spin_box.setMaximumSize(QtCore.QSize(100, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.bundle_n1_spin_box.setFont(font)
        self.bundle_n1_spin_box.setReadOnly(True)
        self.bundle_n1_spin_box.setMinimum(1)
        self.bundle_n1_spin_box.setObjectName("bundle_n1_spin_box")
        self.horizontalLayout_32.addWidget(self.bundle_n1_spin_box)
        self.horizontalLayout_2.addLayout(self.horizontalLayout_32)
        self.horizontalLayout_38 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_38.setObjectName("horizontalLayout_38")
        self.label_21 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label_21.setFont(font)
        self.label_21.setObjectName("label_21")
        self.horizontalLayout_38.addWidget(self.label_21)
        self.bundle_n2_spin_box = QtWidgets.QSpinBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.bundle_n2_spin_box.sizePolicy().hasHeightForWidth())
        self.bundle_n2_spin_box.setSizePolicy(sizePolicy)
        self.bundle_n2_spin_box.setMinimumSize(QtCore.QSize(60, 36))
        self.bundle_n2_spin_box.setMaximumSize(QtCore.QSize(100, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.bundle_n2_spin_box.setFont(font)
        self.bundle_n2_spin_box.setReadOnly(True)
        self.bundle_n2_spin_box.setMinimum(1)
        self.bundle_n2_spin_box.setObjectName("bundle_n2_spin_box")
        self.horizontalLayout_38.addWidget(self.bundle_n2_spin_box)
        self.horizontalLayout_2.addLayout(self.horizontalLayout_38)
        self.horizontalLayout_3.addLayout(self.horizontalLayout_2)
        spacerItem7 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem7)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        self.verticalLayout_3.addLayout(self.verticalLayout_2)
        self.horizontalLayout_4.addLayout(self.verticalLayout_3)
        spacerItem8 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem8)
        self.verticalLayout_4.addLayout(self.horizontalLayout_4)
        self.line_8 = QtWidgets.QFrame(self.centralwidget)
        self.line_8.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_8.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_8.setObjectName("line_8")
        self.verticalLayout_4.addWidget(self.line_8)
        self.verticalLayout_5 = QtWidgets.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.horizontalLayout_31 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_31.setObjectName("horizontalLayout_31")
        self.horizontalLayout_30 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_30.setObjectName("horizontalLayout_30")
        self.horizontalLayout_29 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_29.setObjectName("horizontalLayout_29")
        self.label_17 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(18)
        self.label_17.setFont(font)
        self.label_17.setObjectName("label_17")
        self.horizontalLayout_29.addWidget(self.label_17)
        self.element1_combo_box = QtWidgets.QComboBox(self.centralwidget)
        self.element1_combo_box.setMinimumSize(QtCore.QSize(0, 36))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(18)
        self.element1_combo_box.setFont(font)
        self.element1_combo_box.setObjectName("element1_combo_box")
        self.element1_combo_box.addItem("")
        self.element1_combo_box.addItem("")
        self.element1_combo_box.addItem("")
        self.element1_combo_box.addItem("")
        self.element1_combo_box.addItem("")
        self.horizontalLayout_29.addWidget(self.element1_combo_box)
        self.horizontalLayout_30.addLayout(self.horizontalLayout_29)
        self.horizontalLayout_28 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_28.setObjectName("horizontalLayout_28")
        self.label_18 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(18)
        self.label_18.setFont(font)
        self.label_18.setObjectName("label_18")
        self.horizontalLayout_28.addWidget(self.label_18)
        self.element2_combo_box = QtWidgets.QComboBox(self.centralwidget)
        self.element2_combo_box.setMinimumSize(QtCore.QSize(0, 36))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(18)
        self.element2_combo_box.setFont(font)
        self.element2_combo_box.setObjectName("element2_combo_box")
        self.element2_combo_box.addItem("")
        self.element2_combo_box.addItem("")
        self.element2_combo_box.addItem("")
        self.element2_combo_box.addItem("")
        self.element2_combo_box.addItem("")
        self.horizontalLayout_28.addWidget(self.element2_combo_box)
        self.horizontalLayout_30.addLayout(self.horizontalLayout_28)
        self.horizontalLayout_31.addLayout(self.horizontalLayout_30)
        spacerItem9 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_31.addItem(spacerItem9)
        self.verticalLayout_5.addLayout(self.horizontalLayout_31)
        self.horizontalLayout_24 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_24.setObjectName("horizontalLayout_24")
        self.horizontalLayout_14 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_14.setObjectName("horizontalLayout_14")
        self.elements_bond_label = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(18)
        self.elements_bond_label.setFont(font)
        self.elements_bond_label.setObjectName("elements_bond_label")
        self.horizontalLayout_14.addWidget(self.elements_bond_label)
        self.bond_double_spin_box = QtWidgets.QDoubleSpinBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.bond_double_spin_box.sizePolicy().hasHeightForWidth())
        self.bond_double_spin_box.setSizePolicy(sizePolicy)
        self.bond_double_spin_box.setMinimumSize(QtCore.QSize(150, 36))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(18)
        self.bond_double_spin_box.setFont(font)
        self.bond_double_spin_box.setInputMethodHints(QtCore.Qt.ImhFormattedNumbersOnly)
        self.bond_double_spin_box.setReadOnly(False)
        self.bond_double_spin_box.setDecimals(3)
        self.bond_double_spin_box.setMinimum(0.5)
        self.bond_double_spin_box.setMaximum(10.0)
        self.bond_double_spin_box.setSingleStep(0.01)
        self.bond_double_spin_box.setProperty("value", 1.42)
        self.bond_double_spin_box.setObjectName("bond_double_spin_box")
        self.horizontalLayout_14.addWidget(self.bond_double_spin_box)
        self.horizontalLayout_24.addLayout(self.horizontalLayout_14)
        spacerItem10 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_24.addItem(spacerItem10)
        self.verticalLayout_5.addLayout(self.horizontalLayout_24)
        self.verticalLayout_4.addLayout(self.verticalLayout_5)
        spacerItem11 = QtWidgets.QSpacerItem(20, 8, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_4.addItem(spacerItem11)
        MWNTGenerator.setCentralWidget(self.centralwidget)

        self.retranslateUi(MWNTGenerator)
        self.element1_combo_box.setCurrentIndex(3)
        self.element2_combo_box.setCurrentIndex(3)
        QtCore.QMetaObject.connectSlotsByName(MWNTGenerator)

    def retranslateUi(self, MWNTGenerator):
        _translate = QtCore.QCoreApplication.translate
        MWNTGenerator.setWindowTitle(_translate("MWNTGenerator", "MWNT Generator"))
        self.label_19.setText(_translate("MWNTGenerator", "Generate from..."))
        self.mwnt_Ch_list_radio_button.setText(_translate("MWNTGenerator", "(n, m) List"))
        self.add_Ch_push_button.setText(_translate("MWNTGenerator", "Add (n, m)"))
        self.edit_selected_Ch_push_button.setText(_translate("MWNTGenerator", "Edit Selected"))
        self.remove_selected_Ch_push_button.setText(_translate("MWNTGenerator", "Remove Selected"))
        self.clear_Ch_list_push_button.setText(_translate("MWNTGenerator", "Clear List"))
        self.mwnt_wall_parameters_radio_button.setText(_translate("MWNTGenerator", "Wall Parameters"))
        self.label_24.setText(_translate("MWNTGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">N</span><span style=\" vertical-align:sub;\">walls </span>= </p></body></html>"))
        self.label_25.setText(_translate("MWNTGenerator", "<html><head/><body><p>Wall <span style=\" font-style:italic;\">d</span><span style=\" vertical-align:sub;\">t,min </span>= </p></body></html>"))
        self.min_wall_diameter_double_spin_box.setSuffix(_translate("MWNTGenerator", " Å"))
        self.label_23.setText(_translate("MWNTGenerator", "<html><head/><body><p>Wall <span style=\" font-style:italic;\">d</span><span style=\" vertical-align:sub;\">t,max </span>= </p></body></html>"))
        self.max_wall_diameter_double_spin_box.setSuffix(_translate("MWNTGenerator", " Å"))
        self.label_26.setText(_translate("MWNTGenerator", "<html><head/><body><p>Wall <span style=\" font-family:\'Open Sans,Helvetica Neue,Helvetica,Arial,sans-serif\'; font-size:14px; color:#333333; background-color:#f9f9f9;\">Δ</span><span style=\" font-style:italic;\">d</span><span style=\" vertical-align:sub;\">t,min </span>= </p></body></html>"))
        self.wall_spacing_double_spin_box.setSuffix(_translate("MWNTGenerator", " Å"))
        self.label_5.setText(_translate("MWNTGenerator", "Chiralities"))
        self.random_wall_chiralities_radio_button.setText(_translate("MWNTGenerator", "Random"))
        self.armchair_wall_chiralities_radio_button.setText(_translate("MWNTGenerator", "Armchair"))
        self.zigzag_wall_chiralities_radio_button.setText(_translate("MWNTGenerator", "Zigzag"))
        self.achiral_wall_chiralities_radio_button.setText(_translate("MWNTGenerator", "Achiral"))
        self.chiral_wall_chiralities_radio_button.setText(_translate("MWNTGenerator", "Chiral"))
        self.label_20.setText(_translate("MWNTGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">L</span><span style=\" vertical-align:sub;\"/>= </p></body></html>"))
        self.mwnt_L_double_spin_box.setSuffix(_translate("MWNTGenerator", " Å"))
        self.bundle_generator_check_box.setText(_translate("MWNTGenerator", "Bundle Generator"))
        self.label_8.setText(_translate("MWNTGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">n</span><span style=\" vertical-align:sub;\">1 </span>=</p></body></html>"))
        self.label_21.setText(_translate("MWNTGenerator", "<html><head/><body><p><span style=\" font-style:italic;\">n</span><span style=\" vertical-align:sub;\">2 </span>=</p></body></html>"))
        self.label_17.setText(_translate("MWNTGenerator", "Element 1:"))
        self.element1_combo_box.setItemText(0, _translate("MWNTGenerator", "H"))
        self.element1_combo_box.setItemText(1, _translate("MWNTGenerator", "He"))
        self.element1_combo_box.setItemText(2, _translate("MWNTGenerator", "B"))
        self.element1_combo_box.setItemText(3, _translate("MWNTGenerator", "C"))
        self.element1_combo_box.setItemText(4, _translate("MWNTGenerator", "N"))
        self.label_18.setText(_translate("MWNTGenerator", "Element 2:"))
        self.element2_combo_box.setItemText(0, _translate("MWNTGenerator", "H"))
        self.element2_combo_box.setItemText(1, _translate("MWNTGenerator", "He"))
        self.element2_combo_box.setItemText(2, _translate("MWNTGenerator", "B"))
        self.element2_combo_box.setItemText(3, _translate("MWNTGenerator", "C"))
        self.element2_combo_box.setItemText(4, _translate("MWNTGenerator", "N"))
        self.elements_bond_label.setText(_translate("MWNTGenerator", "<html><head/><body><p>C-C bond =</p></body></html>"))
        self.bond_double_spin_box.setSuffix(_translate("MWNTGenerator", " Å"))

