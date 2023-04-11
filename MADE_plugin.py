# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
import os
#import stat
import sys
import platform
from pathlib import Path

# entry point to PyMOL's API
# pymol.Qt provides the PyQt5 interface, but may support PyQt4
# and/or PySide as well
from pymol.Qt import QtWidgets, QtCore

Plugin_Name_ = "MADE_plugin"
Plugin_Name = "MADE plugin"

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt(Plugin_Name, lambda: run())

dialog = None

def start_plugin():
    # Starts the plugin from the module
    print(("Python version: "+platform.architecture()[0]))
    sys.path.append(os.path.normpath(os.path.join(PLUGIN_DIRECTORY, "module")))
    import MADE_plugin_module
    MADE_plugin_module.main()

class CustomDialog(QtWidgets.QDialog):
    # GUI dialog for plugin install
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setWindowTitle("Choose Directory")
    def find_dir(self):
        folderpath = str(QtWidgets.QFileDialog.getExistingDirectory(CustomDialog(), "Select Directory"))
        if folderpath.endswith("."+Plugin_Name_) or folderpath.endswith("."+Plugin_Name_+os.path.sep):
            pass
        else:
            folderpath=os.path.join(folderpath, "."+Plugin_Name_)
        self.lineEdit.setText(folderpath)

    def set_dir(self):
        global SetDir
        SetDir = self.lineEdit.text()   
        install_dat = open(install_dir, "w")
        install_dat.write(SetDir)
        install_dat.close()
        InstallDialog.close()
            
    def setupUI(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(640, 148)
        Dialog.setMaximumSize(QtCore.QSize(16777215, 210))
        self.gridLayout = QtWidgets.QGridLayout(Dialog)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 2)
        self.lineEdit = QtWidgets.QLineEdit(Dialog)
        self.lineEdit.setObjectName("lineEdit")
        self.gridLayout.addWidget(self.lineEdit, 1, 0, 1, 1)
        HOME_DIRECTORY=os.path.expanduser('~')
        DEFAULT_DIRECTORY=os.path.join(HOME_DIRECTORY, "."+Plugin_Name_)
        self.lineEdit.setText(DEFAULT_DIRECTORY)
        self.pushButton = QtWidgets.QPushButton(Dialog)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 1, 1, 1, 1)
        self.pushButton.clicked.connect(self.find_dir)
        self.label_2 = QtWidgets.QLabel(Dialog)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1)
        self.pushOK = QtWidgets.QPushButton(Dialog)
        self.pushOK.setObjectName("pushOK")
        self.gridLayout.addWidget(self.pushOK, 3, 1, 1, 1)
        self.pushOK.clicked.connect(self.set_dir)
        

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)
        
    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        self.label.setText(_translate("Dialog", "Please choose install directory for {}:".format(Plugin_Name)))
        self.pushButton.setText(_translate("Dialog", "Pick directory"))
        self.label_2.setText(_translate("Dialog", "*Input .{} folder.".format(Plugin_Name_)))
        self.pushOK.setText(_translate("Dialog", "OK"))

def run_custom_dialog():
    global InstallDialog
    InstallDialog=QtWidgets.QDialog()
    InstallDialog.ui = CustomDialog()
    InstallDialog.ui.setupUI(InstallDialog)
    InstallDialog.show()
    InstallDialog.exec_()

def NoInstallSetting():
    # .MADE_plugin_installdir.txt not present, open custom dialog and create it
    global SetDir
    run_custom_dialog()
    if not SetDir=="":
        install_dat = open(install_dir, "w")
        install_dat.write(SetDir)
        install_dat.close()
    else:
        QtWidgets.QMessageBox.about(dialog, Plugin_Name+" Install Warning", "Please choose an install directory!")
        print("Please choose an install directory!")
    
def run():
    global PLUGIN_DIRECTORY, UI_DIRECTORY, install_dir
    # .MADE_plugin_installdir.txt in system home, it points to the install directory
    install_dir=os.path.join(Path.home(), ".{}_installdir.txt".format(Plugin_Name_))

            
    if os.path.isfile(install_dir) == True:
        install_dat = open(install_dir, "r")
        if os.path.isdir(install_dat.readline()):
            install_dat.close()
        else:
            install_dat.close()
            os.remove(install_dir)
        
    try:
        if os.path.isfile(install_dir) == False:
            NoInstallSetting()        
        install_dat = open(install_dir, "r")
        PLUGIN_DIRECTORY=install_dat.readline()
        install_dat.close()
        UI_DIRECTORY=os.path.join(PLUGIN_DIRECTORY,"UI")

    except:
        if os.path.exists(install_dir):
            os.remove(install_dir)
        QtWidgets.QMessageBox.about(dialog, Plugin_Name+" Install Warning", "Something went wrong!\nPlease choose an install directory!")
        print("Something went wrong!\nPlease choose an install directory!")


    print("Initialising {} ...".format(Plugin_Name))
    print(PLUGIN_DIRECTORY)
    global running, status
    status = "start"

    try:
        sys.path.remove('')
    except:
        pass
    
    running = False
    main()

def main():
    # start the plugin
    global status

    if status == "start":
        print("Starting "+Plugin_Name)
        start_plugin()

