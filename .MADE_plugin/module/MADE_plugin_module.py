# Copyright Notice
# ================
#
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
#
# ----------------------------------------------------------------------
# This PyMOL Plugin is Copyright (C) 2022 by Marko Jukic <marko.jukic@um.si>
#
#                        All Rights Reserved
#
# Permission to use, copy and distribute
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------


# -*- coding: utf-8 -*-
# MADE plugin (MAcromolecular DEnsity and Structure Analysis)
# written for python 3.8.x
# Not Yet Described at PyMOL wiki: /
# Author : Marko Jukic
# Date: 2022
# License: UM FKKT Laboratory of Physical Chemistry and Chemical Thermodynamics



#
from __future__ import absolute_import
from __future__ import print_function
import sys, re, gzip, os, urllib.request, urllib.parse, urllib.error, urllib.request, urllib.error, urllib.parse
import stat
import pymol, math, multiprocessing
import subprocess
import numpy as np
from pymol.cgo import *
from pymol import cmd
from pymol import stored
import time
from glob import glob
from ftplib import FTP
from itertools import groupby
from pathlib import Path
import json

from pymol.Qt import QtWidgets, QtGui, QtCore
from pymol.Qt.utils import loadUi

Plugin_Name_ = "MADE_plugin"
Plugin_Name = "MADE plugin"

# SKLEARN----------------------------------------------------------------------=
platform=sys.platform
default_python = "python"
try:
    lib_loc = os.path.dirname(os.path.dirname(os.__file__))
    if platform == "win32":
        default_python_win = os.path.join(lib_loc,"python.exe")
        if os.path.isfile(default_python_win):
            default_python = default_python_win
    else:
        default_python = "python3"
except:
    pass
try:
    from sklearn.cluster import DBSCAN
    print("Sklearn module is installed")
except:
    #from pip._internal import main as pip
    print("Sklearn module has not been installed yet. \nProgram will now install Scikit-learn module")
    time.sleep(1)
    comd = f"{default_python} -m pip install scikit-learn"
    print(comd)
    pro = subprocess.Popen(comd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True, universal_newlines=True)
    stdout,stderr = pro.communicate()
    print(stdout)
    #pymol.cmd.do("pip install scikit-learn")
    #pip(['install', 'scikit-learn'])
    from sklearn.cluster import DBSCAN

#  Timers
time_last = time.time()

def timer_start():
    global time_last
    time_last = time.time()

def time_since_timer_start():
    return round(time.time() - time_last,3)


dict_hetatm_typ_coords_master = {}
dict_hetatm_clusters_master = {}
hetatm_typ_max_consv = {}
hetatm_typ_nr = {}
DBSCAN_eps_of_last_calc = 0

# ---------------------------------------INITIALIZE-----------------------------
install_dir=os.path.join(Path.home(), ".{}_installdir.txt".format(Plugin_Name_))

install_datoteka = open(install_dir, "r")
PLUGIN_DIRECTORY=install_datoteka.readline()
install_datoteka.close()
MODULE_DIRECTORY=os.path.join(PLUGIN_DIRECTORY, "module")
UI_DIRECTORY=os.path.join(PLUGIN_DIRECTORY,"UI")
versionFile = os.path.join(PLUGIN_DIRECTORY, "version.txt")
SETTINGS_DIRECTORY=os.path.join(PLUGIN_DIRECTORY, "settings")
# ------------------------------------------------------------------------------
MHL_dir, report_dir = "",""
def main():
    global default_dl_path, default_dl_site, dl_online_sett, dl_machine_sett, url, check_path, MHL_dir, report_dir, probis_exec, TMalign_exec, gplus_exec, DeepAlign_exec
    print("Settings file for instalation directory location: ", install_dir)
    default_dl_path = os.path.join(PLUGIN_DIRECTORY, Plugin_Name_ + "_DB")
    default_dl_site = "https://cdn.rcsb.org/resources/sequence/clusters"
    dl_online_sett=os.path.join(SETTINGS_DIRECTORY, "DB_Download_From.txt")
    dl_machine_sett=os.path.join(SETTINGS_DIRECTORY, "DB_Download_To.txt")
    dl=open(dl_online_sett, "r") 
    url= dl.read()
    dl.close()
    if url=="":
        dl=open(dl_online_sett, "w")
        dl.write(default_dl_site)
        dl.close()
        dl=open(dl_online_sett, "r") 
        url= dl.read()
        dl.close()
    machine = open (dl_machine_sett, "r")
    check_path = machine.read()
    machine.close()
    if check_path=="":
        machine = open (dl_machine_sett, "w")
        machine.write(default_dl_path)
        machine.close()
        machine = open (dl_machine_sett, "r")
        check_path = machine.read()
        machine.close()

    settings.get_system()

    if platform == "win32":
        probis_exec  = os.path.join(check_path, "probis.exe")
        TMalign_exec = os.path.join(check_path, "TMalign.exe")
        gplus_exec = os.path.join(check_path, "gplus.exe")
        DeepAlign_exec = os.path.join(check_path, "DeepAlign.exe")
    else:
        probis_exec  = os.path.join(check_path, "probis")
        TMalign_exec = os.path.join(check_path, "TMalign")
        gplus_exec = os.path.join(check_path, "gplus")
        DeepAlign_exec = os.path.join(check_path, "DeepAlign")

    GUI().run()


# Running subprocesses
def Run_Subprocess(args):
    args_str = ""
    if type(args) == list:
        for arg in args:
            args_str += str(arg) + " "
    elif type(args) == str:
        args_str = args
    print(args_str)
    pro = subprocess.Popen(args_str,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True, universal_newlines=True)
    stdout,stderr = pro.communicate()
    print(stdout)
    return stdout


################# Proggres Bar, Worker


class PopUpProgressB(QtWidgets.QWidget):
    #Progress bar
    to_run = ""
    def __init__(self):
        super().__init__()
        self.label_2 = QtWidgets.QLabel(self)
        self.label_2.setObjectName("label_2")
        
        self.pbar = QtWidgets.QProgressBar(self)
        self.pbar.setGeometry(30, 40, 500, 75)
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.label_2)
        self.layout.addWidget(self.pbar)
        self.setLayout(self.layout)
        self.setGeometry(300, 300, 550, 100)
        self.setWindowTitle('Progress Bar')

        

        self.worker = Worker()
        self.thread = QtCore.QThread()
        self.worker.prog_sig.connect(self.on_count_changed)
        self.worker.signal_make_msg_box.connect(self.make_msg_box)
        self.worker.prog_msg.connect(self.set_label_prog_msg)
        self.worker.prog_msg.connect(self.set_label_prog_msg)
        self.worker.analysis_GUI.connect(Analysis.Analysis_GUI_interaction)
        self.worker.find_GUI.connect(ClusterComplexManipulation.Find_GUI_interaction)
        self.worker.Filter_HETATM_GUI.connect(HETATM_clustering.Filter_HETATM_GUI_interaction)
        self.worker.Filter_HETATM_GUI_types.connect(HETATM_clustering.Filter_HETATM_GUI_interaction_types)
        self.worker.Download_GUI.connect(BindingSites.get_binding_sites)

        self.worker.moveToThread(self.thread)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.hide) 
        self.thread.started.connect(self.worker.Worker_Run)
        

    def start_progress(self): 
        self.worker.to_run = self.to_run      
        self.thread.start()

    def show_progress(self):
        self.worker.prog_sig.emit(0)
        self.show()

    def on_count_changed(self, value):
        self.pbar.setValue(value)
    
    def set_label_prog_msg(self,prog_msg):
        self.label_2.setText(prog_msg)


    def make_msg_box(self,msg):
        QtWidgets.QMessageBox.about(dialog, msg.split("|")[0],msg.split("|")[1])

def Stop_Worker():
    global popup_bar
    print("Stopping current task")
    popup_bar.worker.exit_flag = True
    popup_bar.worker.finished.emit()
    

class Worker(QtCore.QObject):
    #Worker for long tasks with a progress bar
    finished = QtCore.pyqtSignal()
    prog_sig = QtCore.pyqtSignal(int)
    prog_msg = QtCore.pyqtSignal(str)
    analysis_GUI = QtCore.pyqtSignal(object)
    find_GUI = QtCore.pyqtSignal(str)
    Filter_HETATM_GUI = QtCore.pyqtSignal()
    Filter_HETATM_GUI_types = QtCore.pyqtSignal()
    Download_GUI = QtCore.pyqtSignal()
    to_run = ""
    exit_flag = False

    signal_make_msg_box = QtCore.pyqtSignal(str)
    
    
    def Worker_Run(self):
        self.exit_flag = False
        if self.to_run == "download_complexes":
            self.prog_msg.emit("Downloading Complexes:")
            ClusterComplexManipulation.download_complexes(self)
        elif self.to_run == "analyze_hetams":
            Analysis.analyze_hetams(self)
        elif self.to_run == "Filter_HETATM":
            HETATM_clustering.Filter_HETATM(self)
        elif self.to_run == "get_cluster_complexes":
            self.prog_msg.emit("Finding chains in sequence similarity cluster:")
            ClusterComplexManipulation.get_cluster_complexes(self)
        elif self.to_run == "Setup_DB":
            self.prog_msg.emit("Downloading Sequence identity cluster files")
            Database_setup.Setup_DB(self)






# global reference to avoid garbage collection of our dialog
dialog = None
custom = False
dict_hetatm_typ_coords = {}
analyzesed_complexes_dict = {}

red00 = QtGui.QColor('#ffffff')
red01 = QtGui.QColor('#ffe6e6')
red02 = QtGui.QColor('#ffcccc')
red03 = QtGui.QColor('#ffb3b3')
red04 = QtGui.QColor('#ff9999')
red05 = QtGui.QColor('#ff8080')
red06 = QtGui.QColor('#ff6666')
red07 = QtGui.QColor('#ff4d4d')
red08 = QtGui.QColor('#ff3333')
red09 = QtGui.QColor('#ff1a1a')
red10 = QtGui.QColor('#ff0000')
red_list = [red00,red01,red02,red03,red04,red05,red06,red07,red08,red09,red10]

class GUI:
    def run(self):
        '''
        Open our custom dialog
        '''
        global dialog
    
        dialog = self.make_dialog()
    
        dialog.show()
    
    def make_dialog(self):    
        # create a new Window
        global dialog
        dialog = QtWidgets.QMainWindow()
        dialog.threadpool = QtCore.QThreadPool()
        # populate the Window from our *.ui file which was created with the Qt Designer
        uifile = os.path.join(UI_DIRECTORY, 'MADE_pluginMain.ui')
        self.form = loadUi(uifile, dialog)
        
        # Functions interacting with the GUI
        def dl_changed():
            line_input=dialog.LineDlFrom.text()
            if line_input == url:
                dialog.PushSetDL.setEnabled(False)
            else:
                dialog.PushSetDL.setEnabled(True)
        def machine_changed():
            line_input=dialog.LineDlTo.text()
            #print(check_path)
            if  line_input== check_path:
                dialog.PushSetMachine.setEnabled(False)
            else:
                dialog.PushSetMachine.setEnabled(True)
        def CheckCompare_changed():
            checked = dialog.CheckCompare.isChecked()
        
            if checked == True:
                dialog.CheckWater.setChecked(False)
                dialog.CheckAnalyze.setChecked(False)
                dialog.CheckWater.setEnabled(False)
            else:
                dialog.CheckAnalyze.setChecked(True)

            if checked == True:
                target_complex = dialog.LineProtein.text()
                if len(target_complex) == 4 and os.path.isfile(target_complex+".pdb"):
                    BindingSites.get_binding_sites()
                else:
                    dialog.ListIdentify.clear()

        def CheckAnalyse_changed():
            checked = dialog.CheckAnalyze.isChecked()
            
            if checked == True:
                dialog.CheckCompare.setChecked(False)
                dialog.CheckWater.setEnabled(True)

                target_complex = dialog.LineProtein.text()
                if len(target_complex) == 4 and os.path.isfile(target_complex+".pdb"):
                    BindingSites.get_binding_sites()
                else:
                    dialog.ListIdentify.clear()
            else:
                dialog.CheckWater.setEnabled(False)
                dialog.CheckCompare.setChecked(True)

        def CheckWater_changed():
            target_complex = dialog.LineProtein.text()
            if len(target_complex) == 4 and os.path.isfile(target_complex+".pdb"):
                BindingSites.get_binding_sites()
            else:
                dialog.ListIdentify.clear()

        
        def LineProtein_change(text):
            target_complex = dialog.LineProtein.text()
            if len(target_complex) == 4 and os.path.isfile(target_complex+".pdb"):
                BindingSites.get_binding_sites()
            else:
                dialog.ListIdentify.clear()

        '''
        def CheckPyMolAlign_change():
            if dialog.CheckPyMolAlign.isChecked():
                dialog.CheckPyMolSuper.setChecked(False)

        def CheckPyMolSuper_change():
            if dialog.CheckPyMolSuper.isChecked():
                dialog.CheckPyMolAlign.setChecked(False)
        '''

        def MHL_name_var_change():
            try:
                dialog.LabelPresent.clear()
                curr_item  = dialog.ListIdentify.currentItem().text().split(" ")[0]
                target_complex = dialog.LineProtein.text().lower()
                ComboMMseqs2_list = [30,40,50,70,90,95,100,"Custom"]
                clustering_cutoff = ComboMMseqs2_list[dialog.ComboMMseqs2.currentIndex()]
                superposition_tool = dialog.ComboAlignTool.currentText()

                MHL_name = "master_hetatm_list_{}_{}_{}_{}.txt".format(target_complex,curr_item,clustering_cutoff,superposition_tool)
                dialog.LineMHLName.setText(MHL_name)
            except:
                pass


        def LineMHLName_change(text):
            dialog.CheckUseCurrentMHL.setEnabled(False)
            if os.path.isfile(os.path.join(MHL_dir,dialog.LineMHLName.text())):
                dialog.LabelPresent.setText("Master Heteroatom list already present")
                dialog.CheckUseCurrentMHL.setEnabled(True)
            else:
                dialog.CheckUseCurrentMHL.setChecked(False)

                
        def CheckAllowed_change():
            if dialog.CheckAllowed.isChecked():
                dialog.CheckNotAllowed.setChecked(False) 
                allowed_list = Find.find_allowed_list()

                dialog.LabelAllowedNot.setText("List of allowed HETATMs: ")
                dialog.ListAllowedNot.clear()

                for het in allowed_list:
                    dialog.ListAllowedNot.addItem(het)
            else:
                dialog.LabelAllowedNot.setText("Displaying all HETATMs")
                dialog.ListAllowedNot.clear()


        def CheckNotAllowed_change():            
            if dialog.CheckNotAllowed.isChecked():
                dialog.CheckAllowed.setChecked(False) 
                not_allowed_list = Find.find_not_allowed_list()
                dialog.LabelAllowedNot.setText("List of not allowed HETATMs: ")
                dialog.ListAllowedNot.clear()
                for het in not_allowed_list:
                    dialog.ListAllowedNot.addItem(het)
            else:
                dialog.LabelAllowedNot.setText("Displaying all HETATMs")
                dialog.ListAllowedNot.clear()

        def Tabs_change(ind):
            #curr_tab = dialog.Tabs.currentIndex()
            #print("curr_tab",curr_tab,ind)
            global dict_hetatm_typ_coords_master
            if ind == 1 and dict_hetatm_typ_coords_master:
                HETATM_clustering.run_Filter_HETATM()
                


        def Select_HETATM_type():
            
            HETATM_clustering.Filter_HETATM_GUI_interaction()


        def reset_report_list():
            global report_list_1
            report_list_1 = []


        def protein_change():
            global custom
            custom=False
        
        global popup_bar
        popup_bar = PopUpProgressB()

        # hook up button callbacks
        self.form.PushCustom.clicked.connect(custom_disk_file_get.load_file)
        self.form.PushFind.clicked.connect(ClusterComplexManipulation.run_get_cluster_complexes)
        self.form.PushDownlad.clicked.connect(ClusterComplexManipulation.run_download_complexes)
        self.form.PushDownlad.setEnabled(False)
        #self.form.PushIdentify.clicked.connect(BindingSites.get_binding_sites)
        self.form.PushGo.clicked.connect(Analysis.run_analyze_hetams)
        self.form.PushSetup.clicked.connect(Database_setup.run_Setup_DB)
        self.form.PushDisplay.clicked.connect(pyMOLinterface.pyMOL_display_cluster)
        self.form.PushBSite.clicked.connect(pyMOLinterface.pyMOL_bsite_cluster)
        self.form.PushContacts.clicked.connect(pyMOLinterface.PyMOL_close_resi_contacts)
        self.form.PushChain.clicked.connect(pyMOLinterface.pyMOL_chain_box)
        self.form.PushFetch.clicked.connect(pyMOLinterface.pyMOL_fetch_system)
        self.form.LineDlTo.textChanged.connect(machine_changed)
        self.form.LineDlFrom.textChanged.connect(dl_changed)
        self.form.LineDlFrom.setText(url)
        self.form.LineDlTo.setText(check_path)
        self.form.PushSetDL.clicked.connect(settings.set_download_from)
        self.form.PushSetMachine.clicked.connect(settings.set_download_to)
        self.form.PushCurrentDL.clicked.connect(settings.current_download_from)
        self.form.PushCurrentMachine.clicked.connect(settings.current_download_to)
        self.form.PushDefault.clicked.connect(settings.default)
        self.form.LineProtein.textEdited.connect(protein_change)

        #self.form.PushCluster.clicked.connect(HETATM_clustering.run_Filter_HETATM)
        self.form.CheckAllowed.stateChanged.connect(CheckAllowed_change)
        self.form.CheckNotAllowed.stateChanged.connect(CheckNotAllowed_change)
        self.form.ListIdentify.itemSelectionChanged.connect(MHL_name_var_change)
        self.form.ComboAlignTool.currentIndexChanged.connect(MHL_name_var_change)
        self.form.ComboMMseqs2.currentIndexChanged.connect(MHL_name_var_change)
        self.form.LineMHLName.textChanged.connect(LineMHLName_change)
        self.form.LineProtein.textChanged.connect(LineProtein_change)
        #self.form.CheckPyMolSuper.stateChanged.connect(CheckPyMolSuper_change)
        #self.form.CheckPyMolAlign.stateChanged.connect(CheckPyMolAlign_change)
        self.form.CheckAnalyze.stateChanged.connect(CheckAnalyse_changed)
        self.form.CheckCompare.stateChanged.connect(CheckCompare_changed)
        self.form.CheckWater.stateChanged.connect(CheckWater_changed)


        self.form.PushStop.clicked.connect(Stop_Worker)


        dialog.CheckAllowed.setChecked(True)
        dialog.CheckAllowed.setChecked(False)
        dialog.CheckUseCurrentMHL.setEnabled(False)


        settings.read_settings_from_file()


        self.form.ComboAlignTool.currentIndexChanged.connect(settings.write_settings_to_file)
        self.form.SpinRadius.valueChanged.connect(settings.write_settings_to_file)
        self.form.SpinDB.valueChanged.connect(settings.write_settings_to_file)
        self.form.SpinBsiteRadius.valueChanged.connect(settings.write_settings_to_file)
        self.form.SpinAddBsiteRad.valueChanged.connect(settings.write_settings_to_file)
        self.form.CheckAllowed.stateChanged.connect(settings.write_settings_to_file)
        self.form.CheckNotAllowed.stateChanged.connect(settings.write_settings_to_file)
        #self.form.CheckPyMolAlign.stateChanged.connect(settings.write_settings_to_file)
        #self.form.CheckPyMolSuper.stateChanged.connect(settings.write_settings_to_file)
        
        self.form.LineProtein.textChanged.connect(reset_report_list)
        self.form.ComboAlignTool.currentIndexChanged.connect(reset_report_list)

        self.form.Tabs.tabBarClicked.connect(Tabs_change)
        self.form.ListHetatms.itemSelectionChanged.connect(Select_HETATM_type)


        if os.path.isdir(check_path) == False:
            # First start
            dialog.CheckAligned.setEnabled(False)
            dialog.CheckAnalyze.setEnabled(False)
            dialog.CheckCompare.setEnabled(False)
            dialog.CheckDebye.setEnabled(False)
            dialog.CheckWater.setEnabled(False)
            dialog.ComboMMseqs2.setEnabled(False)
            dialog.LineProtein.setEnabled(False)
            dialog.ListIdentify.setEnabled(False)
            dialog.PushCustom.setEnabled(False)
            dialog.PushFind.setEnabled(False)
            dialog.PushGo.setEnabled(False)
            #dialog.PushIdentify.setEnabled(False)
            dialog.CheckKeep.setEnabled(False)
            dialog.ListCalculated.setEnabled(False)
            dialog.PushBSite.setEnabled(False)
            dialog.PushChain.setEnabled(False)
            dialog.PushContacts.setEnabled(False)
            dialog.PushDisplay.setEnabled(False)
            dialog.PushFetch.setEnabled(False)

            #dialog.PushCluster.setEnabled(False)
            dialog.Tabs.setCurrentIndex(2)

            QtWidgets.QMessageBox.about(dialog, Plugin_Name + " IMPORTANT NOTE", "On the first start of the plugin, please setup the executable/binary DB and database folder. This can be done in 'Settings' tab by clicking on the 'SETUP DATABASE' button! This needs to be done only once.")
            print("Please setup executable/binary DB and database folder first using 'Setup Database' in Settings tab!")
        else:
            global MHL_dir, report_dir
            os.chdir(check_path)
            MHL_dir = os.path.join(check_path,"MHL")
            if not os.path.isdir(MHL_dir):
                os.mkdir(MHL_dir)
            report_dir = os.path.join(check_path,"Reports")
            if not os.path.isdir(report_dir):
                os.mkdir(report_dir)
            Database_setup.file_checks(popup_bar.worker)
        return dialog
    


# ----------------------------------------------------FUNCTIONS-----------------
class settings:
    def set_download_from():
        global url
        set_from = dialog.LineDlFrom.text()
        dl=open(dl_online_sett, "w") 
        dl.write(set_from)
        dl.close()
        dl=open(dl_online_sett, "r") 
        url=dl.read()
        dl.close()        
    def set_download_to():
        global check_path
        set_to = dialog.LineDlTo.text()
        machine = open (dl_machine_sett, "w")
        machine.write(set_to)
        machine.close()
        machine = open (dl_machine_sett, "r")
        check_path = machine.read()
        machine.close()
        line_input=dialog.LineDlTo.text()
        #print(check_path)
        if  line_input== check_path:
            dialog.PushSetMachine.setEnabled(False)
        else:
            dialog.PushSetMachine.setEnabled(True)

        if os.path.isdir(check_path) == False:
            dialog.CheckAligned.setEnabled(False)
            dialog.CheckAnalyze.setEnabled(False)
            dialog.CheckCompare.setEnabled(False)
            dialog.CheckDebye.setEnabled(False)
            dialog.CheckWater.setEnabled(False)
            dialog.ComboMMseqs2.setEnabled(False)
            dialog.LineProtein.setEnabled(False)
            dialog.ListIdentify.setEnabled(False)
            dialog.PushCustom.setEnabled(False)
            dialog.PushFind.setEnabled(False)
            dialog.PushGo.setEnabled(False)
            #dialog.PushIdentify.setEnabled(False)
            dialog.CheckKeep.setEnabled(False)
            dialog.ListCalculated.setEnabled(False)
            dialog.PushBSite.setEnabled(False)
            dialog.PushChain.setEnabled(False)
            dialog.PushContacts.setEnabled(False)
            dialog.PushDisplay.setEnabled(False)
            dialog.PushFetch.setEnabled(False)

            #dialog.PushCluster.setEnabled(False)
            dialog.Tabs.setCurrentIndex(2)
            QtWidgets.QMessageBox.about(dialog, Plugin_Name + " Database Warning", "You have changed the Database Machine directory!\nPlease setup executable/binary DB and database folder first using 'Setup Database' in Settings tab!")
            print("You have changed the Database Machine directory!\nPlease setup executable/binary DB and database folder first using 'Setup Database' in Settings tab!")
        else:
            os.chdir(check_path)
            global probis_exec,TMalign_exec,gplus_exec,DeepAlign_exec,report_dir,MHL_dir
            if platform == "win32":
                probis_exec  = os.path.join(check_path, "probis.exe")
                TMalign_exec = os.path.join(check_path, "TMalign.exe")
                gplus_exec = os.path.join(check_path, "gplus.exe")
                DeepAlign_exec = os.path.join(check_path, "DeepAlign.exe")
            else:
                probis_exec  = os.path.join(check_path, "probis")
                TMalign_exec = os.path.join(check_path, "TMalign")
                gplus_exec = os.path.join(check_path, "gplus")
                DeepAlign_exec = os.path.join(check_path, "DeepAlign")
            MHL_dir = os.path.join(check_path,"MHL")
            if not os.path.isdir(MHL_dir):
                os.mkdir(MHL_dir)
            report_dir = os.path.join(check_path,"Reports")
            if not os.path.isdir(report_dir):
                os.mkdir(report_dir)
            global popup_bar
            Database_setup.file_checks(popup_bar.worker)
    def current_download_from():
        dialog.LineDlFrom.setText(url)
    def current_download_to():
        dialog.LineDlTo.setText(check_path)
    def default():
        dialog.LineDlFrom.setText(default_dl_site)
        dialog.LineDlTo.setText(default_dl_path)
    def get_system():
        global platform
        platform=sys.platform
        
    def write_settings_to_file():
        with open(os.path.join(SETTINGS_DIRECTORY,"Settings.txt"),"w") as set_file:
            #print("Writing Settings:")

            set_file.write("superpos_tool={}\n".format(dialog.ComboAlignTool.currentText()))
            set_file.write("hetatm_rad={}\n".format(dialog.SpinRadius.value() ))
            set_file.write("max_dist={}\n".format(dialog.SpinDB.value() ))
            set_file.write("Local_Bsite_rad={}\n".format(dialog.SpinBsiteRadius.value() ))
            set_file.write("Add_Analysis_rad={}\n".format(dialog.SpinAddBsiteRad.value() ))
            set_file.write("use_allowed={}\n".format(dialog.CheckAllowed.isChecked() ))
            set_file.write("use_not_allowed={}\n".format(dialog.CheckNotAllowed.isChecked() ))
            #set_file.write("PyMol_align={}\n".format(dialog.CheckPyMolAlign.isChecked() ))
            #set_file.write("PyMol_super={}\n".format(dialog.CheckPyMolSuper.isChecked() ))
    
    def read_settings_from_file():
        with open(os.path.join(SETTINGS_DIRECTORY,"Settings.txt"),"r") as set_file:
            lines = set_file.readlines()
            print("Curretnt Settings:\n")

            for line in lines:
                line_split = line.strip().split("=")
                if "superpos_tool" in line:
                    print("superpos_tool",line_split[1])
                    ind = dialog.ComboAlignTool.findText(line_split[1].replace(" ",""))
                    if ind < 0: ind = 0
                    dialog.ComboAlignTool.setCurrentIndex(ind)
                    
                    
                elif "hetatm_rad" in line:
                    print("hetatm_rad",float(line_split[1]))
                    dialog.SpinRadius.setValue(float(line_split[1]))
                elif "max_dist" in line:
                    print("max_dist",float(line_split[1]))
                    dialog.SpinDB.setValue(float(line_split[1]))
                elif "Local_Bsite_rad" in line:
                    print("Bsite_rad",float(line_split[1]))
                    dialog.SpinBsiteRadius.setValue(float(line_split[1]))
                elif "Add_Analysis_rad" in line:
                    print("Add_Analysis_rad",float(line_split[1]))
                    dialog.SpinAddBsiteRad.setValue(float(line_split[1]))
                elif "use_allowed" in line:
                    bl = False
                    if "true" in line_split[1].lower(): 
                        bl = True
                    print("use_allowed",bl)
                    dialog.CheckAllowed.setChecked(bl)
                elif "use_not_allowed" in line:
                    bl = False
                    if "true" in line_split[1].lower(): 
                        bl = True
                    print("use_not_allowed",bl)
                    dialog.CheckNotAllowed.setChecked(bl)
                '''elif "PyMol_align" in line:
                    bl = False
                    if "true" in line_split[1].lower(): 
                        bl = True
                    print("PyMol_align",bl)
                    dialog.CheckPyMolAlign.setChecked(bl)

                elif "PyMol_super" in line:
                    bl = False
                    if "true" in line_split[1].lower(): 
                        bl = True
                    print("PyMol_super",bl)
                    dialog.CheckPyMolSuper.setChecked(bl)
                '''

        


                    
    

class custom_disk_file_get:
    # custom PDB file!
    @staticmethod
    def load_file():
        global custom,filename,check_path
        filename, _filter = QtWidgets.QFileDialog.getOpenFileName(dialog, "Open..." + "PDB File", ".", "PDB File (*.pdb)")
        if filename:
            try:
                print("File read.")
                print((str(filename)))
                
                cp_cmd = "cp"
                if platform == "win32":
                    cp_cmd = "copy"
                
                dest_file = os.path.join(check_path,Path(filename).stem+".pdb").replace("/","\\")
                file_to_cp = str(filename).replace("/","\\")
                Run_Subprocess("{} {} {}".format(cp_cmd,file_to_cp,dest_file))

            except:
                QtWidgets.QMessageBox.about(dialog, Plugin_Name + " Warning", "File error or \nFile not Found! \nPlease investigate!")
            dialog.LineProtein.setText(Path(filename).stem)
            
            custom=True


class Database_setup:
    """PDB database setup"""

    list_cluster_files = ["clusters-by-entity-30.txt", "clusters-by-entity-40.txt", "clusters-by-entity-50.txt", "clusters-by-entity-70.txt","clusters-by-entity-90.txt", "clusters-by-entity-95.txt", "clusters-by-entity-100.txt"]


    @staticmethod
    def run_Setup_DB():
        global popup_bar
        popup_bar.to_run = "Setup_DB"
        
        popup_bar.show_progress()
        popup_bar.start_progress()


    @staticmethod
    def Setup_DB(worker):

        current_path = str(os.getcwd())
        if current_path == check_path:
            pass
        else:
            try:
                os.mkdir(check_path)
            except:
                pass
            os.chdir(check_path)


        """download the cluster database"""
        print("Setting up the cluster database!")
        prog_count = 0
        for cluster_file in Database_setup.list_cluster_files:
            print("Fetching",cluster_file)
            url_file = url+"/"+cluster_file
            urllib.request.urlretrieve(url_file, cluster_file)
            prog_count += 1
            prog = int(prog_count/len(Database_setup.list_cluster_files)*100)
            worker.prog_sig.emit(prog)
            if worker.exit_flag:
                break
        
        print(("Database setup finished" + "\t thank you to the RCSB protein data bank"))
        superpos_methods = ["probis","TMalign","DeepAlign","gplus"]
        
        
        worker.prog_msg.emit("Dowloading Superposition method executables")

        prog_count = 0
        worker.prog_sig.emit(0)
        for sp_m in superpos_methods:
            if worker.exit_flag:
                break
            Database_setup.fetch_superpos(sp_m)
            prog_count += 1
            prog = int(prog_count/len(superpos_methods)*100)
            worker.prog_sig.emit(prog)

        global probis_exec,TMalign_exec,gplus_exec,DeepAlign_exec,report_dir,MHL_dir

        if platform == "linux" or platform == "win32":
            pass
        else:

            worker.signal_make_msg_box.emit(Plugin_Name+" Warning"+"|"+"System error:\nOperating on a non-compatible operating system!")

        if platform == "win32":
            probis_exec  = os.path.join(check_path, "probis.exe")
            TMalign_exec = os.path.join(check_path, "TMalign.exe")
            gplus_exec = os.path.join(check_path, "gplus.exe")
            DeepAlign_exec = os.path.join(check_path, "DeepAlign.exe")
        else:
            probis_exec  = os.path.join(check_path, "probis")
            TMalign_exec = os.path.join(check_path, "TMalign")
            gplus_exec = os.path.join(check_path, "gplus")
            DeepAlign_exec = os.path.join(check_path, "DeepAlign")
        MHL_dir = os.path.join(check_path,"MHL")
        if not os.path.isdir(MHL_dir):
            os.mkdir(MHL_dir)
        report_dir = os.path.join(check_path,"Reports")
        if not os.path.isdir(report_dir):
            os.mkdir(report_dir)

        if not os.path.isfile("clusters_custom.txt"):
            open("clusters_custom.txt", 'a').close()

        dialog.CheckAligned.setEnabled(True)
        dialog.CheckAnalyze.setEnabled(True)
        dialog.CheckCompare.setEnabled(True)
        dialog.CheckDebye.setEnabled(True)
        dialog.CheckWater.setEnabled(False)
        dialog.ComboMMseqs2.setEnabled(True)
        dialog.LineProtein.setEnabled(True)
        dialog.ListIdentify.setEnabled(True)
        dialog.PushCustom.setEnabled(True)
        dialog.PushFind.setEnabled(True)
        dialog.PushGo.setEnabled(True)
        #dialog.PushIdentify.setEnabled(True)
        dialog.CheckKeep.setEnabled(True)
        dialog.ListCalculated.setEnabled(True)
        dialog.PushBSite.setEnabled(True)
        dialog.PushChain.setEnabled(True)
        dialog.PushContacts.setEnabled(True)
        dialog.PushDisplay.setEnabled(True)
        dialog.PushFetch.setEnabled(True)

        #dialog.PushCluster.setEnabled(True)
        Database_setup.file_checks(worker)
        worker.finished.emit()
        return None


    @staticmethod
    def fetch_superpos(superpos):
        #downloading superposition method executables from our repo
        print("Fetching {} executable".format(superpos))
        if platform == "win32":
            superpos += ".exe"
            if superpos == "probis.exe":
                urllib.request.urlretrieve("https://gitlab.com/Vid-Ra/test/-/raw/files/libgsl.dll", "libgsl.dll")
                urllib.request.urlretrieve("https://gitlab.com/Vid-Ra/test/-/raw/files/libgslcblas.dll", "libgslcblas.dll")
            elif superpos == "gplus.exe":
                urllib.request.urlretrieve("https://gitlab.com/Vid-Ra/test/-/raw/files/libgcc_s_dw2-1.dll", "libgcc_s_dw2-1.dll")
                urllib.request.urlretrieve("https://gitlab.com/Vid-Ra/test/-/raw/files/libstdc++-6.dll", "libstdc++-6.dll")

        
        urllib.request.urlretrieve("https://gitlab.com/Vid-Ra/test/-/raw/files/{}".format(superpos), superpos) # temp URL ##################################################
        if platform == "linux":
            if os.access(superpos, os.X_OK) == True:
                print("{} executable ok".format(superpos))
            else:
                st = os.stat(superpos)
                os.chmod(superpos, st.st_mode | stat.S_IEXEC)
                print("{} executable permission set".format(superpos))

    @staticmethod
    def file_checks(worker):
        """Checking database files"""
        global probis_exec,TMalign_exec,DeepAlign_exec,gplus_exec
        print("\n" +Plugin_Name+": checking database!")
        check = ""
        all_ok = True
        for file in Database_setup.list_cluster_files:
            print((file + " present " + str(os.path.isfile(file))))
            if str(os.path.isfile(file)) == "False":
                print(Plugin_Name+": please setup DB of rscb cluster files")
                
                worker.signal_make_msg_box.emit(Plugin_Name+" Warning"+"|"+ "please setup DB of rscb cluster files using the SETUP DATABASE button in Settings tab on plugin GUI")
                dialog.Tabs.setCurrentIndex(2)
                all_ok = False
                break
            else:
                check += "a"
        superpos_methods = ["probis","TMalign","DeepAlign","gplus"]
        superpos_methods_execs = [probis_exec,TMalign_exec,DeepAlign_exec,gplus_exec]
        
        if check == "aaaaaaaaaaa":
            print(Plugin_Name+ ": Database OK")
        
        popup_bool = False
        for i,superpos_method in enumerate(superpos_methods):
            if not os.path.isfile(superpos_methods_execs[i]):
                popup_bool = True
                print(Plugin_Name+": please ensure {} is downloaded correctly using the SETUP DATABASE button in Settings tab on plugin GUI".format(superpos_method))
        if popup_bool:
            all_ok = False
            dialog.Tabs.setCurrentIndex(2)
            worker.signal_make_msg_box.emit(Plugin_Name+" Warning"+"|"+"Please ensure the superposition method executables are downloaded correctly using the SETUP DATABASE button in Settings tab on plugin GUI")

        if platform == "linux":
            for i,superpos_method in enumerate(superpos_methods):
                try:
                    if os.access(superpos_methods_execs[i], os.X_OK) == True:
                        print("{} executable ok".format(superpos_method))
                    else:
                        st = os.stat(superpos_methods_execs[i])
                        os.chmod(superpos_methods_execs[i], st.st_mode | stat.S_IEXEC)
                        print("{} executable permission set".format(superpos_method))
                except:
                    all_ok = False
        if all_ok:  
            dialog.Tabs.setCurrentIndex(0)
            dialog.CheckAligned.setEnabled(True)
            dialog.CheckAnalyze.setEnabled(True)
            dialog.CheckCompare.setEnabled(True)
            dialog.CheckDebye.setEnabled(True)
            dialog.CheckWater.setEnabled(False)
            dialog.ComboMMseqs2.setEnabled(True)
            dialog.LineProtein.setEnabled(True)
            dialog.ListIdentify.setEnabled(True)
            dialog.PushCustom.setEnabled(True)
            dialog.PushFind.setEnabled(True)
            dialog.PushGo.setEnabled(True)
            #dialog.PushIdentify.setEnabled(True)
            dialog.CheckKeep.setEnabled(True)
            dialog.ListCalculated.setEnabled(True)
            dialog.PushBSite.setEnabled(True)
            dialog.PushChain.setEnabled(True)
            dialog.PushContacts.setEnabled(True)
            dialog.PushDisplay.setEnabled(True)
            dialog.PushFetch.setEnabled(True)

            #dialog.PushCluster.setEnabled(True)



# ------------------------------------------------------------------------------


analyzesed_complexes_dict_temp = {}
class ThreadWorker_get_chain_list(QtCore.QRunnable):
    # worker class for multithreaded querying the pdb for the chains belonging to an entity
    def __init__(self,comp,entity):
        super().__init__()
        self.comp = comp
        self.entity = entity

    def run(self):
        global analyzesed_complexes_dict_temp,fin_thread_worker_counter,popup_bar,thread_worker_counter

        if popup_bar.worker.exit_flag:
            return
        chain_ls = ClusterComplexManipulation.get_chain_list_from_pdb_entity_id(self.comp,self.entity)
        if self.comp not in analyzesed_complexes_dict_temp:
            analyzesed_complexes_dict_temp[self.comp] = chain_ls
        else:
            analyzesed_complexes_dict_temp[self.comp] += chain_ls
        fin_thread_worker_counter += 1
        prog = int(fin_thread_worker_counter/thread_worker_counter*100)
        popup_bar.worker.prog_sig.emit(prog)

class ClusterComplexManipulation:
    """download and identification of cluster complexes from RCSB"""

    @staticmethod
    def Find_GUI_interaction(Find_GUI):
        dialog.LineFind.setText(Find_GUI)
        dialog.PushDownlad.setEnabled(True)

    @staticmethod
    def find_analyzesed_complexes_dict(name,filename,worker):
        # searches filename (clusters-by-entity-xx.txt or clusters_custom.txt)
        # finds the analyzesed_complexes_dict, dictionary of analyzesed_complexes_dict[PDB_ID] = [list_of_similar_chains]

        if dialog.threadpool.maxThreadCount() > 16:
            dialog.threadpool.setMaxThreadCount(16)
        
        vzorec_name =re.compile(name,re.IGNORECASE) 
        lin_split = []
        
        try:
            clus_str = "MMseqs2"
            if filename == "clusters_custom.txt":
                clus_str = "Custom"
            else:
                seq_cutoff = filename.replace("clusters-by-entity-","").replace(".txt","")
                clus_str = seq_cutoff + "% sequence identity MMseqs2"


            with open(filename , "rt") as infile:
                for lin in infile.readlines():
                    match = vzorec_name.search(lin)
                    
                    if match:
                        lin_split_temp = lin.strip().split(' ')
                        for comp_chain in lin_split_temp:
                            match2 = vzorec_name.search(comp_chain)
                            
                            if match2:
                                
                                try:
                                    comp,entity = comp_chain.split("_")
                                except:
                                    comp = ""
                                    spl = comp_chain.split("_")
                                    for i in range(len(spl)-1):
                                        comp += spl[i]
                                        if i != len(spl)-2:
                                            comp += "_"
                                    entity = spl[-1]
                                if len(comp) == 4:
                                    lin_split = lin_split_temp
                                    break

        except:
            print(filename,"does not exist!")
            return
        
        global fin_thread_worker_counter,thread_worker_counter
        fin_thread_worker_counter = 0
        thread_worker_counter = 0
        target_complex = dialog.LineProtein.text()

        worker.find_GUI.emit(str(len(lin_split)) + " complexes in {} cluster of {}.".format(clus_str,target_complex))


        global analyzesed_complexes_dict_temp
        analyzesed_complexes_dict_temp = {}
        for comp_chain in lin_split:
            try:
                comp,entity = comp_chain.split("_")
            except:
                comp = ""
                spl = comp_chain.split("_")
                for i in range(len(spl)-1):
                    comp += spl[i]
                    if i != len(spl)-2:
                        comp += "_"
                entity = spl[-1]

            if worker.exit_flag:
                
                break

            comp = comp.lower()
            entity = entity.upper()
            if len(comp )== 4:
                thread_worker_counter += 1
                thread_worker = ThreadWorker_get_chain_list(comp,entity)
                dialog.threadpool.start(thread_worker)

        
        while fin_thread_worker_counter != thread_worker_counter:
            if worker.exit_flag:
                break
           
        return analyzesed_complexes_dict_temp
    @staticmethod
    def get_chain_list_from_pdb_entity_id(pdb_id,entity_id):
        # queries the rcsb pdb for chains belonging to the entity
        try:
            url = "https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}".format(pdb_id.upper(),entity_id)
            response = urllib.request.urlopen( url ) 
            js = json.load(response)
            chain_ls = js["entity_poly"]["pdbx_strand_id"].split(",")
        except:
            print("Error with getting chains for {} {}".format(pdb_id.upper(),entity_id))
            print(sys.exc_info())
            chain_ls = []

        return chain_ls


    @staticmethod
    def get_cluster_unique_list(target_complex, sele_name, seq_id,worker):

        target_complex = dialog.LineProtein.text()
        print("target_complex",target_complex)
        try:
        
            global analyzesed_complexes_dict
            if sele_name == "clusters":
                filename = sele_name + seq_id + ".txt"
                analyzesed_complexes_dict =  ClusterComplexManipulation.find_analyzesed_complexes_dict(target_complex,filename,worker)
                clus_str = "Custom"

                #worker.find_GUI.emit(str(len(analyzesed_complexes_dict)) + " complexes in cluster")
                print("Using custom cluster of complexes, reading data form clusters_custom.txt in the database directory")
                print("Found " + str(len(analyzesed_complexes_dict)) +  " of complexes in cluster")

            else:
                filename = sele_name + seq_id + ".txt"
                print("filename",filename)
                analyzesed_complexes_dict =  ClusterComplexManipulation.find_analyzesed_complexes_dict(target_complex,filename,worker)
                clus_str = seq_id + "% sequence identity MMseqs2"

                print("\nComplexes and chains: ")
                print(analyzesed_complexes_dict)

                print("using precalculated sqeuence identity clusters from the PDB with MMseqs2 clustering:")
                print("""Steinegger, M., Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35, 1026–1028 (2017). https://doi.org/10.1038/nbt.3988""")
                print("Found " + str(len(analyzesed_complexes_dict)) +  " of complexes in cluster")

            n_chains = sum([len(chains) for i,chains in analyzesed_complexes_dict.items()])
            worker.find_GUI.emit(str(len(analyzesed_complexes_dict)) + " complexes in {} cluster of {}. {} homologous chains.".format(clus_str,target_complex,n_chains))
            worker.finished.emit()

        
        except:
            worker.signal_make_msg_box.emit(Plugin_Name+" Warning |Invalid PDB ID or \nDatabase File not Found! \nPlease investigate!")
            print(sys.exc_info())
            worker.finished.emit()

    @staticmethod
    def run_get_cluster_complexes():
        
        global popup_bar
        popup_bar.to_run = "get_cluster_complexes"
        
        popup_bar.show_progress()
        popup_bar.start_progress()

    @staticmethod
    def get_cluster_complexes(worker):

        target_complex = dialog.LineProtein.text()

        if len(target_complex) < 4 and not custom:
            worker.signal_make_msg_box.emit(Plugin_Name+" Warning |Invalid PDB ID!")
            worker.finished.emit()
            
        elif len(target_complex) > 4 and not custom:
            worker.signal_make_msg_box.emit(Plugin_Name+" Warning |Invalid PDB ID!")
            worker.finished.emit()
        else:

            if dialog.ComboMMseqs2.currentText() == "Custom Cluster":
                selected_sequence = "_custom"
                ClusterComplexManipulation.get_cluster_unique_list(target_complex, "clusters", selected_sequence,worker)

            else: 
                selected_sequence = dialog.ComboMMseqs2.currentText().replace("%","")
                ClusterComplexManipulation.get_cluster_unique_list(target_complex, "clusters-by-entity-", selected_sequence,worker)

        not_present = []
        for cmpl in analyzesed_complexes_dict.keys():
            if not os.path.isfile(str(cmpl).lower() + ".pdb"):
                not_present.append(cmpl)
        
        if not_present:
            print("{} Complexes not present, press the Download Complexes button.".format(len(not_present)))
            print("List of missing complexes:",not_present)
        else:
            print("All complexes in the seqeuence identity cluster already downloaded.")
                

    @staticmethod
    def run_download_complexes():
        global popup_bar
        popup_bar.to_run = "download_complexes"
        
        popup_bar.show_progress()
        popup_bar.start_progress()


    @staticmethod
    def Get_PDB_file(pdb_id,ftp_wwpdb_server):
        # downloads a pdb file 
        try:
            dir = '/pub/pdb/data/structures/divided/pdb/' + str(pdb_id[1:3]).lower() + "/"
            #ftp_wwpdb_server.cwd()

            zipped_file = "pdb" + str(pdb_id).lower() + ".ent.gz"
            file_uncompressed = str(pdb_id).lower() + ".pdb"
        
            local_file = open(zipped_file, 'wb')
            ftp_wwpdb_server.retrbinary('RETR ' +dir+ zipped_file, local_file.write)
            local_file.close()
            # unzip

            with gzip.open(zipped_file, 'rb') as compressed_f:
                local_file_uncompressed = open(file_uncompressed, 'wb')
                local_file_uncompressed.write(compressed_f.read())
                local_file_uncompressed.close()
            os.remove(zipped_file)
            return 0
        except:

            print(sys.exc_info())
            print("Problem Downlading PDB ID",pdb_id)
            return 1
    @staticmethod
    def download_complexes(worker):


        global fin_thread_worker_counter,thread_worker_counter,ftp_wwpdb_server
        fin_thread_worker_counter,thread_worker_counter = 0,0
        
        try:
            
                # -----

            print("Downloading files... \nThis may take a while depending on the number of complexes")
            ftp_wwpdb_server = FTP("ftp.wwpdb.org")
            ftp_wwpdb_server.login()
            print(len(analyzesed_complexes_dict),analyzesed_complexes_dict)
            for complex in analyzesed_complexes_dict.keys():


                if not os.path.isfile(str(complex).lower() + ".pdb"):
                    ret = ClusterComplexManipulation.Get_PDB_file(complex,ftp_wwpdb_server)
                    if not ret:
                        print(("Downloaded complex " +str(complex).lower() +", "+ str(thread_worker_counter+1) + " out of " + str(len(analyzesed_complexes_dict))))
                else:

                    print("Complex Existing "+str(complex).lower() +", "+ str(thread_worker_counter+1) + " out of " + str(len(analyzesed_complexes_dict)))
                fin_thread_worker_counter += 1

                prog = int(fin_thread_worker_counter/len(analyzesed_complexes_dict)*100)
                worker.prog_sig.emit(prog)
                thread_worker_counter += 1
                if worker.exit_flag:
                    
                    ftp_wwpdb_server.quit()
                    break

            ftp_wwpdb_server.quit()
            print("Download of complexes finished")
            worker.signal_make_msg_box.emit(Plugin_Name+ " Message" + "|Download of complexes finished!")
            
            
            worker.Download_GUI.emit()

            worker.finished.emit()
            return None
        except:
            worker.signal_make_msg_box.emit(Plugin_Name+" Warning"+"|"+ "Check all settings, \nPlease investigate!")
        worker.finished.emit()
        


class BindingSites:
    """definiraj binding site"""

    bsite_unique_centers = []

    @staticmethod
    def get_binding_sites():
        # finds binding sites in a pdb file
        try:
            dialog.ListIdentify.clear()
    
            water_sel = dialog.CheckWater.isChecked()
            chain_sel = dialog.CheckCompare.isChecked()
            if custom == True:
                target_complex_2 = filename #dialog.LineProtein.text()
            else:
                target_complex_2 = (dialog.LineProtein.text().lower() + ".pdb")
            
            vzorec_2 = re.compile("^" + "ATOM")
            vzorec_3 = re.compile("^" + "HETATM")
            vzorec_4 = re.compile("HOH")
            # LISTS---------------------------------------------------------
            list_heteratoms = []
            list_atoms = []
            list_waters = []
            list_binding_sites = []
            list_water_binding_sites = []
            list_chains = []
            list_chains_konc = []
            warning_lista = []
            bsite_unique = []
            # Results in lists --------------------------------
            try:
                with open(target_complex_2, "rt") as infile:
                    for linenumber, line in enumerate(infile):
                        if vzorec_3.search(line) != None:
                            if vzorec_4.search(line) != None:
                                list_waters.append(line.rstrip('\n'))
                            else:
                                list_heteratoms.append(line.rstrip('\n'))
                        elif vzorec_2.search(line) != None:
                            list_atoms.append(line.rstrip('\n'))
            except OSError:
                QtWidgets.QMessageBox.about(dialog, Plugin_Name + " Warning", "File not found \n Please Investigate!")
                print("File not found!, please investigate.")
                return
            
            global list_atoms_xyzchain
            list_atoms_xyzchain = []
    
            for entry in list_atoms:
                test = []
                try:
                    # ATOM    477  CG2 ILE A  78       6.540   0.762  34.941  1.00 13.42           C
                    entry[30:38]
                    test.append(entry[30:38])
                    test.append(entry[38:46])
                    test.append(entry[46:54])
                    test.append(str(entry[21]))
                    list_atoms_xyzchain.append(test)
                except:
                    pass
    
    
            if water_sel == False:
                for linija in list_heteratoms:
                    unique_binding_site = linija[17:20].strip(" ")+"."+linija[22:26].strip(" ")+"."+linija[21].strip(" ")
                    list_binding_sites.append(unique_binding_site)
                    bsite_unique.append([unique_binding_site, linija[30:38].strip(" "), linija[38:46].strip(" "), linija[46:54].strip(" ")])
            else:
                for linija in list_heteratoms+list_waters:
                    unique_binding_site = linija[17:20].strip(" ")+"."+linija[22:26].strip(" ")+"."+linija[21].strip(" ")
                    list_binding_sites.append(unique_binding_site)
                    bsite_unique.append([unique_binding_site, linija[30:38].strip(" "), linija[38:46].strip(" "), linija[46:54].strip(" ")])
    
    
    
            # Binding site prep ---------------------------------------------------
            for key, group in groupby(bsite_unique, lambda x: x[0]):
                bsx = []
                bsy = []
                bsz = []
                for el in group:
                    bsx.append(float(el[1]))
                    bsy.append(float(el[2]))
                    bsz.append(float(el[3]))
    
                # name of bsite, axerage x, average y, average z, min x, max x, min y, max y, min z, max z
                BindingSites.bsite_unique_centers.append([key, sum(bsx)/len(bsx), sum(bsy)/len(bsy), sum(bsz)/len(bsz), min(bsx), max(bsx), min(bsy), max(bsy), min(bsz), max(bsz)])
    
    
    
            # waters -----------------------------------------------------
            for linija in list_waters:
                water_binding_sites = linija[17:20].strip(" ")+"."+linija[22:26].strip(" ")+"."+linija[21].strip(" ")
                list_water_binding_sites.append(water_binding_sites)
    
            # other atoms ---------------------------------------------------
            for linija in list_atoms:
                atom_site = []
                try:
                    atom_site.append(str(linija[21]))
                    atom_site.append(int(linija[22:26]))
                    list_chains.append(atom_site)
                except:
                    pass
    
            # whole chain ---------------------------------------------------
            for key, group in groupby(list_chains, lambda x: x[0]):
                temp_group = []
                #residue_number = 0
    
                for el in group:
                    # counting unique residues
                    if el not in temp_group:
                        temp_group.append(el)
                    else:
                        pass
                    ins_str = str(temp_group[-1][0]) + " chain with " + str(len(temp_group)) + " residues"
    
    
                list_chains_konc.append(ins_str)

            #-------------------------------------------------------------------
    
            if chain_sel == False:
                for chainelement in list_chains_konc:
                    if float(chainelement.split()[3]) < 30.0:
                        temp_str = "ONLY " + str(chainelement.split()[3]) + " residues found in chain " + str(chainelement.split()[0])
                        warning_lista.append(temp_str)
                    else:
                        pass
    
                if len(warning_lista) >= 1:
                    print( "Try to compare whole chains instead of individual binding sites at the short chain location.")
                else:
                    pass
    
                if water_sel == False:
                    for entry in sorted(list(set(list_binding_sites))):
                        dialog.ListIdentify.addItem(entry)
                else:
                    for entry in sorted(list(set(list_binding_sites))):
                        dialog.ListIdentify.addItem(entry)
                    for entry in sorted(list(set(list_water_binding_sites))):
                        dialog.ListIdentify.addItem(entry)
    
            else:
                for entry in list_chains_konc:
                    dialog.ListIdentify.addItem(entry)
    
            if dialog.ListIdentify.count() == 0:
                if dialog.CheckCompare.isChecked() == False:
                    QtWidgets.QMessageBox.about(dialog, Plugin_Name +" Message",
                                                "No binding sites found for analysis.\nPlease use Compare Whole Chain option or define Water as binding site")
                    print("No binding sites found for analysis.\nPlease use Compare Whole Chain option or define Water as binding site")
                else:
                    QtWidgets.QMessageBox.about(dialog, Plugin_Name+" Message",
                                                "No chains found to compare.")
                    print("No chains found to compare.")
    
            return None
        # flow out except-------------------------------------------------------
        except:
            QtWidgets.QMessageBox.about(dialog, Plugin_Name + " Warning", "Invalid PDB ID or \nDatabase File not Found! \n\nPlease investigate!")
            print((sys.exc_info()))


# report lista 1----------------------------------------------------------------
report_list_1 = []
report_list_2 = []

# ------------------------------------------------------------------------------



#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#
#                                                                                            Alignment functions
#
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################

# GNU Terry Pratchett
Probis_cite = '''Performing protein superposition using ProBiS: http://insilab.org/probis-algorithm/
Konc, J., Miller, B. T., Stular, T., Lesnik, S., Woodcock, H. L., Brooks, B. R., and Janezic, D. ProBiS-CHARMMing: Web Interface for Prediction and Optimization of Ligands in Protein Binding Sites. J. Chem. Inf. Model., 2015, 55, 2308-2314.
'''
gplus_cite = '''Performing protein superposition using GANGSTA+: https://github.com/guerler/gplus
Guerler, A. and Knapp, E.-W. (2008), Novel protein folds and their nonsequential structural analogs. Protein Science, 17: 1374-1382. https://doi.org/10.1110/ps.035469.108
'''

Deep_align_cite = '''Performing protein superposition using DeepAlign: https://github.com/realbigws/DeepAlign
Sheng Wang, Jianzhu Ma, Jian Peng and Jinbo Xu. PROTEIN STRUCTURE ALIGNMENT BEYOND SPATIAL PROXIMITY Scientific Reports, 3, 1448, (2013)
'''

Pymol_align_cite = '''Performing protein superposition using the align command in PyMOL: https://pymolwiki.org/index.php/Align'''
Pymol_super_cite = '''Performing protein superposition using the super command in PyMOL: https://pymolwiki.org/index.php/Super'''

TM_align_cite = '''Performing protein superposition using TM-align: https://zhanggroup.org/TM-align/
Y. Zhang, J. Skolnick, TM-align: A protein structure alignment algorithm based on TM-score, Nucleic Acids Research, 33: 2302-2309 (2005)
'''

class Align_methods:
    # methods for protein alignment

    @staticmethod
    def return_appropriate_chain_list(aligned_complex,target_complex,chain_list,one_or_multiple):
        # filters target complex and returns only one chain if single aligned chain per PDB entry checked
        chain_list_out = chain_list
        if not os.path.isfile(aligned_complex+".pdb"):
            print(aligned_complex+".pdb not present")
            return []
        if aligned_complex == target_complex:
            return []
        try:
            if one_or_multiple:
                chain_list_out = [chain_list[0]]
        except:
            pass
        return chain_list_out



    @staticmethod
    def PyMol_save_chain_het(aligned_complex,chain,filename):
        # saves a specific chain and all the heteratoms regardless of chain 
        cmd.load(aligned_complex+".pdb")
        cmd.select("chain_hetatm","({} and chain {}) or ({} and  hetatm)".format(aligned_complex,chain,aligned_complex))
        cmd.save(filename,"chain_hetatm")
        cmd.delete("chain_hetatm")
        cmd.delete(aligned_complex)




    @staticmethod
    def align_DeepAlign_one_pdb(target_complex,bsite_selection,whole_chain_compare_selection,aligned_complex,chain_id,chain_sel_bool,filename_target):
        # Running DeepAlign for one pdb
        global DeepAlign_exec
        filename_align  = aligned_complex+"_{}.chain.pdb".format(chain_id)
        cmd.load(aligned_complex+".pdb")
        cmd.select("chain_hetatm","({} and chain {})".format(aligned_complex,chain_id))
        cmd.save(filename_align,"chain_hetatm")
        if not chain_sel_bool:
            filnam = target_complex+bsite_selection.split(".")[2]+"_"+aligned_complex+chain_id+".0.rota.pdb"
        else:
            filnam = target_complex+whole_chain_compare_selection+"_"+aligned_complex+chain_id+".0.rota.pdb"

        ''' if platform == "win32":
            subprocess.run(args = [DeepAlign_exec,filename_align,filename_target,"-o","Align_Deep"])

        else:
            subprocess.run(args = [DeepAlign_exec,filename_align,filename_target,"-o","Align_Deep"])'''
        Run_Subprocess(args = [DeepAlign_exec,filename_align,filename_target,"-o","Align_Deep"])

        aligned_complex_sele = "rota_{}_{}".format(aligned_complex,chain_id)
        filename_het_chain  = aligned_complex+"_{}_het.chain.pdb".format(chain_id)
        cmd.select("chain_hetatm","({} and chain {}) or ({} and  hetatm)".format(aligned_complex,chain_id,aligned_complex))
        cmd.save(filename_het_chain,"chain_hetatm")
        cmd.delete("chain_hetatm")
        cmd.delete(aligned_complex)
        cmd.load(filename_het_chain,aligned_complex_sele)


        try:
                    
            # reading the rotation matrx
            rot_score_filename = "Align_Deep.score"

            rot_mat = []
            with open(rot_score_filename,"r") as stdout_file:
                lines = stdout_file.readlines()

                
                for ln in range(4,7):
                    mat_ln = lines[ln].strip().split(" ")
                    while "" in mat_ln:
                        mat_ln.remove("")

                    rot_mat.append(float(mat_ln[0]))

                    for i in range(1,len(mat_ln)):
                        rot_mat.append(float(mat_ln[i]))

            rot_mat.append(0)
            rot_mat.append(0)
            rot_mat.append(0)
            rot_mat.append(0)

            # manual rotaton with PyMOL
            cmd.transform_selection(aligned_complex_sele,rot_mat)
            cmd.save(filnam,aligned_complex_sele)
            cmd.disable(aligned_complex_sele)   

            
        except:
            print("problem with DeepAlign",aligned_complex,chain_id)
            print(sys.exc_info())




    def align_DeepAlign(processors_available_local,bsite_selection,target_complex,chain_selection,analyzesed_complexes_dict,chain_sel_bool,whole_chain_compare_selection,one_or_multiple,worker):
        # alignment with DeepAlign
        cmd.delete("all")
        cmd.load(target_complex+".pdb")
        if not chain_sel_bool:
            Bsite_rad = dialog.SpinBsiteRadius.value()
            cmd.select("target_bsite","byres {} and resi {} and chain {} around {}".format(target_complex,bsite_selection.split(".")[1],bsite_selection.split(".")[2],Bsite_rad)) 

        else:
            cmd.select("target_bsite","{} and chain {}".format(target_complex,whole_chain_compare_selection))
        
        filename_target  = target_complex+"_{}.chain.pdb".format("bsite")
        cmd.save(filename_target,"target_bsite")
        prog_count = 0
        for aligned_complex,chain_list in analyzesed_complexes_dict.items():
            print("Working on complex {}".format(aligned_complex))
            appropriate_chain_list = Align_methods.return_appropriate_chain_list(aligned_complex,target_complex,chain_list,one_or_multiple)
            print("    Chains in sequence identity cluster: {}\n    Analysed chains: {}".format(chain_list,appropriate_chain_list))

            for chain_id in appropriate_chain_list:
                Align_methods.align_DeepAlign_one_pdb(target_complex,bsite_selection,whole_chain_compare_selection,aligned_complex,chain_id,chain_sel_bool,filename_target)
            prog_count += 1
            prog = int(prog_count/len(analyzesed_complexes_dict.keys())*100)
            worker.prog_sig.emit(prog)
            if worker.exit_flag:
                
                break




    @staticmethod
    def align_gplus_one_pdb(target_complex,bsite_selection,whole_chain_compare_selection,aligned_complex,chain_id,chain_sel_bool,filename_target):
        # alignment with GANGSTA+
        global gplus_exec
        filename_align  = aligned_complex+"_{}.chain.pdb".format(chain_id)
        cmd.load(aligned_complex+".pdb")
        cmd.select("chain_hetatm","({} and chain {})".format(aligned_complex,chain_id))
        cmd.save(filename_align,"chain_hetatm")

        if not chain_sel_bool:
            filnam = target_complex+bsite_selection.split(".")[2]+"_"+aligned_complex+chain_id+".0.rota.pdb"
        else:
            filnam = target_complex+whole_chain_compare_selection+"_"+aligned_complex+chain_id+".0.rota.pdb"
        
        stdout_lines = Run_Subprocess("{} {} {}".format(gplus_exec,filename_align,filename_target)).split("\n")
        aligned_complex_sele = "rota_{}_{}".format(aligned_complex,chain_id)
        filename_het_chain  = aligned_complex+"_{}_het.chain.pdb".format(chain_id)

        cmd.select("chain_hetatm","({} and chain {}) or ({} and  hetatm)".format(aligned_complex,chain_id,aligned_complex))
        cmd.save(filename_het_chain,"chain_hetatm")
        cmd.delete("chain_hetatm")
        cmd.delete(aligned_complex)
        cmd.load(filename_het_chain,aligned_complex_sele)

        
        
        try:
            rot_mat = []

                
            for ln in range(14,17):
                mat_ln = stdout_lines[ln].split(" ")
                while "" in mat_ln:
                    mat_ln.remove("")
                mat_ln.remove(mat_ln[0])
                for i in range(1,len(mat_ln)):
                    rot_mat.append(float(mat_ln[i]))
                rot_mat.append(float(mat_ln[0]))

            rot_mat.append(0)
            rot_mat.append(0)
            rot_mat.append(0)
            rot_mat.append(0)


            # manual rotation with PyMOL
            cmd.transform_selection(aligned_complex_sele,rot_mat)
            cmd.save(filnam,aligned_complex_sele)
            cmd.disable(aligned_complex_sele)   
        

        except:
            print("problem with gplus",aligned_complex,chain_id)
            print(sys.exc_info())


    @staticmethod
    def align_gplus(processors_available_local,bsite_selection,target_complex,chain_selection,analyzesed_complexes_dict,chain_sel_bool,whole_chain_compare_selection,one_or_multiple,worker):
        cmd.delete("all")
        cmd.load(target_complex+".pdb")
        if not chain_sel_bool:
            Bsite_rad = dialog.SpinBsiteRadius.value()
            cmd.select("target_bsite","byres (({} and resi {} and chain {} around {}) and polymer.protein)".format(target_complex,bsite_selection.split(".")[1],bsite_selection.split(".")[2],Bsite_rad)) 

        else:
            cmd.select("target_bsite","{} and chain {}".format(target_complex,whole_chain_compare_selection))
        
        filename_target  = target_complex+"_{}.chain.pdb".format("bsite")
        cmd.save(filename_target,"target_bsite")
        prog_count = 0
        for aligned_complex,chain_list in analyzesed_complexes_dict.items():
            print("Working on complex {}".format(aligned_complex))
            appropriate_chain_list = Align_methods.return_appropriate_chain_list(aligned_complex,target_complex,chain_list,one_or_multiple)
            print("    Chains in sequence identity cluster: {}\n    Analysed chains: {}".format(chain_list,appropriate_chain_list))

            for chain_id in appropriate_chain_list:
                Align_methods.align_gplus_one_pdb(target_complex,bsite_selection,whole_chain_compare_selection,aligned_complex,chain_id,chain_sel_bool,filename_target)
            prog_count += 1
            prog = int(prog_count/len(analyzesed_complexes_dict.keys())*100)
            worker.prog_sig.emit(prog)
            if worker.exit_flag:
                
                break


    @staticmethod
    def align_TMalign_one_pdb(target_complex,bsite_selection,whole_chain_compare_selection,aligned_complex,chain_id,chain_sel_bool,filename_target):
        # Aligning one file wiht TM-align
        global TMalign_exec
        filename_align  = aligned_complex+"_{}.chain.pdb".format(chain_id)
        Align_methods.PyMol_save_chain_het(aligned_complex,chain_id,filename_align)
        if not chain_sel_bool:
            filnam = target_complex+bsite_selection.split(".")[2]+"_"+aligned_complex+chain_id+".0.rota.pdb"
        else:
            filnam = target_complex+whole_chain_compare_selection+"_"+aligned_complex+chain_id+".0.rota.pdb"
        
        # run with name Align_TM, output aligned output we are interested is then Align_TM_all_atm_lig.pdb
        if platform == "win32":
            Run_Subprocess(args = [TMalign_exec,filename_align,filename_target,"-o","Align_TM"])

            Run_Subprocess("move Align_TM_all_atm_lig.pdb "+filnam)


        else:
            Run_Subprocess(args = [TMalign_exec,filename_align,filename_target,"-o","Align_TM"])
            Run_Subprocess("mv Align_TM_all_atm_lig.pdb "+filnam)

        aligned_complex_sele = "rota_{}_{}".format(aligned_complex,chain_id)
        try:
            cmd.load(filnam,aligned_complex_sele)
            cmd.disable(aligned_complex_sele)
        except:
            pass


    @staticmethod
    def align_TMalign(processors_available_local,bsite_selection,target_complex,chain_selection,analyzesed_complexes_dict,chain_sel_bool,whole_chain_compare_selection,one_or_multiple,worker):
        cmd.delete("all")
        cmd.load(target_complex+".pdb")
        if not chain_sel_bool:
            Bsite_rad = dialog.SpinBsiteRadius.value()
            cmd.select("target_bsite","byres {} and resi {} and chain {} around {}".format(target_complex,bsite_selection.split(".")[1],bsite_selection.split(".")[2],Bsite_rad)) 

        else:
            cmd.select("target_bsite","{} and chain {}".format(target_complex,whole_chain_compare_selection))
        
        filename_target  = target_complex+"_{}.chain.pdb".format("bsite")
        cmd.save(filename_target,"target_bsite")
        prog_count = 0
        for aligned_complex,chain_list in analyzesed_complexes_dict.items():
            print("Working on complex {}".format(aligned_complex))
            appropriate_chain_list = Align_methods.return_appropriate_chain_list(aligned_complex,target_complex,chain_list,one_or_multiple)
            print("    Chains in sequence identity cluster: {}\n    Analysed chains: {}".format(chain_list,appropriate_chain_list))

            for chain_id in appropriate_chain_list:
                Align_methods.align_TMalign_one_pdb(target_complex,bsite_selection,whole_chain_compare_selection,aligned_complex,chain_id,chain_sel_bool,filename_target)
            prog_count += 1
            prog = int(prog_count/len(analyzesed_complexes_dict.keys())*100)
            worker.prog_sig.emit(prog)
            if worker.exit_flag:
                
                break
    @staticmethod
    def align_Pymol_one_pdb(target_complex,bsite_selection,whole_chain_compare_selection,aligned_complex,chain_id,chain_sel_bool,superposition_tool):
        # aligning one PDB with PyMOL
        aligned_complex_sele = "rota_{}_{}".format(aligned_complex,chain_id)
        cmd.load(aligned_complex+".pdb",aligned_complex_sele)


        if not chain_sel_bool:
            filnam = target_complex+bsite_selection.split(".")[2]+"_"+aligned_complex+chain_id+".0.rota.pdb"
        else:
            filnam = target_complex+whole_chain_compare_selection+"_"+aligned_complex+chain_id+".0.rota.pdb"



        cmd.select("align_bsite","{} and chain {}".format(aligned_complex_sele,chain_id))
        if not chain_sel_bool:
            #if dialog.CheckPyMolAlign.isChecked():
            if "align" in superposition_tool:
                cmd.align("align_bsite","target_bsite")
            else:
                cmd.super("align_bsite","target_bsite") 

        else:
            #if dialog.CheckPyMolAlign.isChecked():
            if "align" in superposition_tool:
                cmd.align("align_bsite","target_bsite")
            else:
                cmd.super("align_bsite","target_bsite") 

        cmd.save(filnam,aligned_complex_sele)

        cmd.disable(aligned_complex_sele)

        cmd.delete("align_bsite")



    @staticmethod
    def align_Pymol(processors_available_local,bsite_selection,target_complex,chain_selection,analyzesed_complexes_dict,chain_sel_bool,whole_chain_compare_selection,one_or_multiple,worker,superposition_tool):
        cmd.delete("all")
        cmd.load(target_complex+".pdb")
        if not chain_sel_bool:
            Bsite_rad = dialog.SpinBsiteRadius.value()
            cmd.select("target_bsite","byres {} and resi {} and chain {} around {}".format(target_complex,bsite_selection.split(".")[1],bsite_selection.split(".")[2],Bsite_rad))

        else:
            cmd.select("target_bsite","{} and chain {}".format(target_complex,whole_chain_compare_selection))

        prog_count = 0
        for aligned_complex,chain_list in analyzesed_complexes_dict.items():
            print("Working on complex {}".format(aligned_complex))
            appropriate_chain_list = Align_methods.return_appropriate_chain_list(aligned_complex,target_complex,chain_list,one_or_multiple)
            print("    Chains in sequence identity cluster: {}\n    Analysed chains: {}".format(chain_list,appropriate_chain_list))

            for chain_id in appropriate_chain_list:
                Align_methods.align_Pymol_one_pdb(target_complex,bsite_selection,whole_chain_compare_selection,aligned_complex,chain_id,chain_sel_bool,superposition_tool)
            prog_count += 1
            prog = int(prog_count/len(analyzesed_complexes_dict.keys())*100)
            worker.prog_sig.emit(prog)
            if worker.exit_flag:
                
                break


    @staticmethod
    def align_ProBiS(processors_available_local,bsite_selection,target_complex_2,chain_selection,analyzesed_complexes_dict,chain_sel,whole_chain_compare_selection,one_or_multiple,worker):
        # Alignment with ProBiS
        worker.prog_msg.emit("Running Superposition part 1")
        global probis_exec
        if platform == "win32": 
            if chain_sel == False:
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-extract", "-bsite", bsite_selection, "-dist 3.0", "-f1",
                        "{}.pdb".format(target_complex_2), "-c1", chain_selection, 
                        " -srffile", "{}.srf".format(target_complex_2)])
            else:
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-extract", "-f1", "{}.pdb".format(target_complex_2), "-c1", 
                        whole_chain_compare_selection, "-srffile", 
                        "{}.srf".format(target_complex_2)])
        else:
            if chain_sel == False:
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-extract", "-bsite", bsite_selection, "-dist 3.0", "-f1",
                        "{}.pdb".format(target_complex_2), "-c1", chain_selection, 
                        " -srffile", "{}.srf".format(target_complex_2)])
            else:
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-extract", "-f1", "{}.pdb".format(target_complex_2), "-c1", 
                        whole_chain_compare_selection, "-srffile", 
                        "{}.srf".format(target_complex_2)])
                


        open("./srfs.txt", 'w').close()
        prot_list = []
        prog_count = 0
        for pdb_id,chain_list in analyzesed_complexes_dict.items():
            print("Working on complex {}".format(pdb_id))
            appropriate_chain_list = Align_methods.return_appropriate_chain_list(pdb_id,target_complex_2,chain_list,one_or_multiple)
            print("    Chains in sequence identity cluster: {}\n    Analysed chains: {}".format(chain_list,appropriate_chain_list))
            for chain_id in appropriate_chain_list:

                if platform == "win32":
                    Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                            "-extract", "-f1", "{}.pdb".format(pdb_id), "-c1", 
                            chain_id, "-srffile", "{}{}.srf".format(pdb_id, chain_id)])
                else:
                    Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                            "-extract", "-f1", "{}.pdb".format(pdb_id), "-c1", 
                            chain_id, "-srffile", "{}{}.srf".format(pdb_id, chain_id)])

                srf_file = open("./srfs.txt", 'a')
                srf_file.write(pdb_id + chain_id + ".srf " + chain_id + "\n")
                srf_file.close()
                prot_list.append(pdb_id + " " + chain_id)
                pass
            prog_count += 1
            prog = int(prog_count/len(analyzesed_complexes_dict.keys())*100)
            worker.prog_sig.emit(prog)
            if worker.exit_flag:
                
                break


        ##############################
        if platform == "win32":
            Run_Subprocess("DEL *.rota.pdb /S")
            Run_Subprocess("DEL AAA_NOSQL.nosql /S")
        else:
            Run_Subprocess("rm ./*.rota.pdb")
            Run_Subprocess("rm ./AAA_NOSQL.nosql")
        ##############################
        
        if chain_sel == False:
            if platform == "win32":
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-surfdb", "-local", "-sfile", "srfs.txt", "-f1", 
                        "{}.srf".format(target_complex_2), "-c1", 
                        chain_selection, "-nosql", "AAA_NOSQL.nosql"])
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-results", "-f1", "{}.pdb".format(target_complex_2), "-c1", 
                        chain_selection, "-nosql", "AAA_NOSQL.nosql", "-json", "AAA_NOSQL.json"])
            else:
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-surfdb", "-local", "-sfile", "srfs.txt", "-f1", 
                        "{}.srf".format(target_complex_2), "-c1", 
                        chain_selection, "-nosql", "AAA_NOSQL.nosql"])
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-results", "-f1", "{}.pdb".format(target_complex_2), "-c1", 
                        chain_selection, "-nosql", "AAA_NOSQL.nosql", "-json", "AAA_NOSQL.json"])
        else:
            if platform == "win32":
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-surfdb", "-sfile", "srfs.txt", "-f1", "{}.srf".format(target_complex_2), 
                        "-c1", whole_chain_compare_selection, "-nosql", "AAA_NOSQL.nosql"])
                
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-results", "-f1", "{}.pdb".format(target_complex_2), "-c1", 
                        whole_chain_compare_selection, "-nosql", "AAA_NOSQL.nosql", "-json", "AAA_NOSQL.json"])
            else:
                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-surfdb", "-sfile", "srfs.txt", "-f1", "{}.srf".format(target_complex_2), 
                        "-c1", whole_chain_compare_selection, "-nosql", "AAA_NOSQL.nosql"])

                Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                        "-results", "-f1", "{}.pdb".format(target_complex_2), "-c1", 
                        whole_chain_compare_selection, "-nosql", "AAA_NOSQL.nosql", "-json", "AAA_NOSQL.json"])
        prog_count = 0
        worker.prog_sig.emit(0)
        worker.prog_msg.emit("Running Superposition part 2")
        for element in prot_list:
            ele0=element.split()[0]
            ele1=element.split()[1]
            
            if chain_sel == False:
                if platform == "win32":
                    Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                            "-align", "-bkeep", "-alno", "0", "-nosql", "AAA_NOSQL.nosql", "-f1", 
                            "{}.pdb".format(target_complex_2), "-c1", chain_selection, 
                            "-f2", "{}.pdb".format(ele0), "-c2", 
                            ele1])
                else:
                    Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                            "-align", "-bkeep", "-alno", "0", "-nosql", "AAA_NOSQL.nosql", "-f1", 
                            "{}.pdb".format(target_complex_2), "-c1", chain_selection, 
                            "-f2", "{}.pdb".format(ele0), "-c2", 
                            ele1])
            else:
                if platform == "win32":
                    Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                            "-align", "-bkeep", "-alno", "0", "-nosql", "AAA_NOSQL.nosql", "-f1", 
                            "{}.pdb".format(target_complex_2), "-c1", 
                            whole_chain_compare_selection, "-f2", 
                            "{}.pdb".format(ele0), "-c2", ele1]) 
                else:
                    Run_Subprocess(args=[probis_exec, "-ncpu", processors_available_local, 
                            "-align", "-bkeep", "-alno", "0", "-nosql", "AAA_NOSQL.nosql", "-f1", 
                            "{}.pdb".format(target_complex_2), "-c1", 
                            whole_chain_compare_selection, "-f2", 
                            "{}.pdb".format(ele0), "-c2", ele1]) 
            prog_count += 1
            prog = int(prog_count/len(prot_list)*100)
            worker.prog_sig.emit(prog)
            if worker.exit_flag:
                
                break












#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#
#                                                                                            finding functions
#
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################

dic_all_hetatm_types_names = {}
dic_all_hetatm_types_formulas = {}
class Find:


    @staticmethod
    def find_hetnames(filename_pdb):
        #finds the names and formulas of heteroatoms form pdb files
        global dic_all_hetatm_types_names
        global dic_all_hetatm_types_formulas

        dic_all_hetatm_types_names_local = {}
        dic_all_hetatm_types_formulas_local = {}
        try:
            fil_pdb = open(filename_pdb,'r')
            lines_pdb = fil_pdb.readlines()
            fil_pdb.close()
        except:
            lines_pdb = []
        vzorec_end_model = re.compile("^ATOM")

        vzorec_hetnam =re.compile("^HETNAM",re.IGNORECASE) 
        vzorec_formul =re.compile("^FORMUL",re.IGNORECASE) 

        for line in lines_pdb:
            if vzorec_end_model.match(line):
                break

            if vzorec_hetnam.match(line):
                typ = line[11:14].strip()
                name = line[15:].strip()

                if typ not in dic_all_hetatm_types_names_local:
                    dic_all_hetatm_types_names_local[typ] = name
                else:
                    dic_all_hetatm_types_names_local[typ] += name

            elif vzorec_formul.match(line):
                typ = line[12:15].strip()
                formul = line[16:].strip()
                if typ not in dic_all_hetatm_types_formulas_local:
                    dic_all_hetatm_types_formulas_local[typ] = formul
                else:
                    dic_all_hetatm_types_formulas_local[typ] +=  formul

        try:
            for typ in dic_all_hetatm_types_names_local.keys():
                if typ not in dic_all_hetatm_types_names:
                    dic_all_hetatm_types_names[typ] = dic_all_hetatm_types_names_local[typ] 
                    dic_all_hetatm_types_formulas[typ] = dic_all_hetatm_types_formulas_local[typ] 
        except:
            pass
    
    @staticmethod
    def find_ion_list():
        # returns the list of heteratoms defined as ions, "List_Ions.txt" in settings directory
        ion_list = []
        with open(os.path.join(SETTINGS_DIRECTORY, "List_Ions.txt"),"r") as fil_ions:
            ion_list = fil_ions.readlines()
        for i in range(len(ion_list)):
            ion_list[i] = ion_list[i].upper().strip()
        return ion_list
    
    @staticmethod
    def find_allowed_list():
        # returns the list of allowed heteroatoms, "List_allowed_hetams.txt" in settings directory
        allowed_list = []
        with open(os.path.join(SETTINGS_DIRECTORY, "List_allowed_hetams.txt"),"r") as fil_allowed:
            allowed_list = fil_allowed.readlines()
        for i in range(len(allowed_list)):
            allowed_list[i] = allowed_list[i].upper().strip()
        return allowed_list


    @staticmethod
    def find_not_allowed_list():
        # returns the list of not allowed heteroatoms, "List_not_allowed_hetams.txt" in settings directory
        not_allowed_list = []
        with open(os.path.join(SETTINGS_DIRECTORY, "List_not_allowed_hetams.txt"),"r") as fil_not_allowed:
            not_allowed_list = fil_not_allowed.readlines()
        for i in range(len(not_allowed_list)):
            not_allowed_list[i] = not_allowed_list[i].upper().strip()
        return not_allowed_list

    @staticmethod
    def find_hetams_from_master_list():
        # returns dict_hetatm_typ_coords from a master heteroatom list file
        dict_hetatm_typ_coords = {}
        bsite_space_check = dialog.CheckAnalyze.isChecked()
        chain_sel = dialog.CheckCompare.isChecked()


        with open(os.path.join(MHL_dir,dialog.LineMHLName.text()), "r") as mhl_file:
            lines = mhl_file.readlines()
            
            for line in lines:
                if line[0] == "#":
                    continue
                x_hetm,y_hetm,z_hetm,B_hetm,serial_nr,chain_iden,hetatm_typ,atom_name,seq_nr,PDB_id = line.strip().split(" ")
                
                
                hetatm_iden = hetatm_typ + "-" + atom_name
                temp_list = [float(x_hetm),float(y_hetm),float(z_hetm)]
                
                if bsite_space_check == True and chain_sel == False:
                    #if SELECTED_SITE[4] - 4 <= temp_list[0] <= SELECTED_SITE[5] + 4:
                        #if SELECTED_SITE[6] - 4 <= temp_list[1] <= SELECTED_SITE[7] + 4:
                            #if SELECTED_SITE[8] - 4 <= temp_list[2] <= SELECTED_SITE[9] + 4:
                    if atom_min_x <= temp_list[0] <= atom_max_x:
                        if atom_min_y <= temp_list[1] <= atom_max_y:
                            if atom_min_z <= temp_list[2] <= atom_max_z:
                                if hetatm_iden not in dict_hetatm_typ_coords:
                                    dict_hetatm_typ_coords[hetatm_iden]  = [[float(x_hetm),float(y_hetm),float(z_hetm)]]
                                else:
                                    dict_hetatm_typ_coords[hetatm_iden].append([float(x_hetm),float(y_hetm),float(z_hetm)])

                if bsite_space_check == False and chain_sel == False:
                    if atom_min_x <= temp_list[0] <= atom_max_x:
                        if atom_min_y <= temp_list[1] <= atom_max_y:
                            if atom_min_z <= temp_list[2] <= atom_max_z:
                                if hetatm_iden not in dict_hetatm_typ_coords:
                                    dict_hetatm_typ_coords[hetatm_iden]  = [[float(x_hetm),float(y_hetm),float(z_hetm)]]
                                else:
                                    dict_hetatm_typ_coords[hetatm_iden].append([float(x_hetm),float(y_hetm),float(z_hetm)])

                if chain_sel == True:
                    if atom_min_x <= temp_list[0] <= atom_max_x:
                        if atom_min_y <= temp_list[1] <= atom_max_y:
                            if atom_min_z <= temp_list[2] <= atom_max_z:
                                if hetatm_iden not in dict_hetatm_typ_coords:
                                    dict_hetatm_typ_coords[hetatm_iden]  = [[float(x_hetm),float(y_hetm),float(z_hetm)]]
                                else:
                                    dict_hetatm_typ_coords[hetatm_iden].append([float(x_hetm),float(y_hetm),float(z_hetm)])
            

        return dict_hetatm_typ_coords 














#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#
#                                                                                            Clustering functions
#
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################


class HETATM_clustering:

    @staticmethod
    def reduce_dict_hetatm_typ_coords_master(print_lists = False):
        global dict_hetatm_typ_coords_master
        # reduces the master heteroatom dictionary with applied allowed/not allowed lists
        dict_hetatm_typ_coords = {}
        bool_list_allowed = dialog.CheckAllowed.isChecked()
        bool_list_not_allowed = dialog.CheckNotAllowed.isChecked()

        allowed_list = []
        not_allowed_list = []
        if bool_list_allowed:
            allowed_list =  Find.find_allowed_list()
            if print_lists:
                print("Allowed list used: ",allowed_list)
        
        if bool_list_not_allowed:
            not_allowed_list = Find.find_not_allowed_list()
            if print_lists:

                print("Not Allowed list used: ",not_allowed_list)

        for hetatm_iden in dict_hetatm_typ_coords_master.keys():
            hetatm_typ, atom_name = hetatm_iden.split("-")
            if bool_list_allowed:
                if hetatm_typ.upper() not in allowed_list:
                    continue
            elif bool_list_not_allowed:
                if hetatm_typ.upper() in not_allowed_list:
                    continue
            dict_hetatm_typ_coords[hetatm_iden] = dict_hetatm_typ_coords_master[hetatm_iden]
        return dict_hetatm_typ_coords


    @staticmethod
    def add_to_dict_hetatm_typ_coords_master(x_hetm,y_hetm,z_hetm,B_hetm,serial_nr,chain_iden,hetatm_typ,atom_name,seq_nr,PDB_id):
        # write an entry to dict_hetatm_typ_coords_master while seaching through superposed files
        global dict_hetatm_typ_coords_master
        bsite_space_check = dialog.CheckAnalyze.isChecked()
        chain_sel = dialog.CheckCompare.isChecked()
        temp_list = [float(x_hetm),float(y_hetm),float(z_hetm)]
        
        hetatm_iden = hetatm_typ + "-" + atom_name

        # checks only for heteroatoms near the binding site/chain 
        if bsite_space_check == True and chain_sel == False:
            #if SELECTED_SITE[4] - 4 <= temp_list[0] <= SELECTED_SITE[5] + 4:
                #if SELECTED_SITE[6] - 4 <= temp_list[1] <= SELECTED_SITE[7] + 4:
                    #if SELECTED_SITE[8] - 4 <= temp_list[2] <= SELECTED_SITE[9] + 4:
            if atom_min_x <= temp_list[0] <= atom_max_x:
                if atom_min_y <= temp_list[1] <= atom_max_y:
                    if atom_min_z <= temp_list[2] <= atom_max_z:                    
                        if hetatm_iden not in dict_hetatm_typ_coords_master:
                            dict_hetatm_typ_coords_master[hetatm_iden]  = [[float(x_hetm),float(y_hetm),float(z_hetm)]]
                        else:
                            dict_hetatm_typ_coords_master[hetatm_iden].append([float(x_hetm),float(y_hetm),float(z_hetm)])

        if bsite_space_check == False and chain_sel == False:
            if atom_min_x <= temp_list[0] <= atom_max_x:
                if atom_min_y <= temp_list[1] <= atom_max_y:
                    if atom_min_z <= temp_list[2] <= atom_max_z:
                        if hetatm_iden not in dict_hetatm_typ_coords_master:
                            dict_hetatm_typ_coords_master[hetatm_iden]  = [[float(x_hetm),float(y_hetm),float(z_hetm)]]
                        else:
                            dict_hetatm_typ_coords_master[hetatm_iden].append([float(x_hetm),float(y_hetm),float(z_hetm)])

        if chain_sel == True:
            if atom_min_x <= temp_list[0] <= atom_max_x:
                if atom_min_y <= temp_list[1] <= atom_max_y:
                    if atom_min_z <= temp_list[2] <= atom_max_z:
                        if hetatm_iden not in dict_hetatm_typ_coords_master:
                            dict_hetatm_typ_coords_master[hetatm_iden]  = [[float(x_hetm),float(y_hetm),float(z_hetm)]]
                        else:
                            dict_hetatm_typ_coords_master[hetatm_iden].append([float(x_hetm),float(y_hetm),float(z_hetm)])



    @staticmethod
    def make_hetatm_iden_clusters(mlist, hetatm_typ):
        # finds all the differend clusters of for hetatm_typ heteroatoms with coordinats in mlist, adds entry to dict_hetatm_clusters_master
        global dict_hetatm_clusters_master
        global entities

        def calculate_clusters_ions(lista, sample_size):
            list_np_array = np.array(lista)
            selected_eps=dialog.SpinDB.value()
            labels3D = DBSCAN(eps=selected_eps, min_samples=sample_size).fit_predict(list_np_array)
            # Number of clusters in labels, ignoring noise if present.
            n_clusters_3D = len(set(labels3D)) - (1 if -1 in labels3D else 0)
            return int(n_clusters_3D)
        
        start_population = 2
        clus_num = calculate_clusters_ions(mlist, start_population)
        cluster_collate_list = []
        max_population = 1
        temp_report = []
        temp_report.append("IDENTIFIED {} CLUSTERS: ".format(hetatm_typ))

        while clus_num >= 1:

            consv = round((float(start_population)/float(entities)), 2)
            # limita ena je overloaded ker je lahko teoreticno prisotnih vec molekul vode na istem mestu v istem kristalu
            # glede na to da je smisel tega orodja v eksperimentalnih podatkih bi bile taksne vode na lokaciji manjsi od 1 A
            # nesmiselne in korigirane s strani kristalografa
            # zato lahko komot v skrajno nenavadno ali eksp-nekorigiranem primeru vrednost consv presega 1
            # taksne primere tukaj reduciramo na vrednost 1 kar pomeni, da je voda na tej lokaciji nastopa v vseh eksperimentalnih entitetah
            if consv > 1:
                consv = 1.0
            else:
                pass

            # report list------------------------------------------------
            text_consv = int(round(consv*10)) * "*"
            st = 10 - len(text_consv)
            text_consv += st * " "
            if str(clus_num) == "1":
                text = (str(clus_num) + " cluster with " + str(start_population)
                        + " {} HETATMS, ".format(hetatm_typ) +  "conservation " + str(consv))
            else:
                text = (str(clus_num) + " clusters with at least " + str(start_population)
                        + " {} HETATMS, ".format(hetatm_typ) +  "conservation " + str(consv))
            num_spaces = 50 - len(text)
            temp_report.append(text +  + num_spaces * " " + "[" + text_consv + "]")
            cluster_collate_list.append([clus_num, start_population, text])
            start_population += 1
            clus_num = calculate_clusters_ions(mlist, start_population)
            max_population = start_population

        dict_hetatm_clusters_master[hetatm_typ] = cluster_collate_list
        return max_population-1,round((float(max_population-1)/float(entities)), 2)



    @staticmethod
    def display_cluster_info_hetatm_typ(hetatm_typ):
        # display all the clusters of hetatm_typ in the calculated clusters display
        global dict_hetatm_clusters_master
        global entities
        
        cluster_collate_list = dict_hetatm_clusters_master[hetatm_typ]
        temp_collate = []
        d=0
        counter = 0
        for cluster_num, start_population, list_text in reversed(cluster_collate_list):
            d=d+1
            
            
            if cluster_num not in temp_collate:
                if counter == 0:
                    list_text = list_text.replace("at least ","")
                dialog.ListCalculated.addItem(list_text)
                temp_collate.append(cluster_num)
                counter += 1
            else:
                pass

        calculated_items = []
        for x in range(dialog.ListCalculated.count()):
            calculated_items.append(dialog.ListCalculated.item(x))
        # red colour according to cluster conservation
        for i, listbox_entry in enumerate(calculated_items):
            try:
                consv = float(listbox_entry.text().split()[7])
            except:
                consv = float(listbox_entry.text().split()[9])
            if 0.9 <= consv :
                dialog.ListCalculated.item(i).setBackground(red10)
            elif 0.8 <=consv < 0.9:
                dialog.ListCalculated.item(i).setBackground(red09)
            elif 0.7 <=consv < 0.8:
                dialog.ListCalculated.item(i).setBackground(red08)
            elif 0.6 <=consv < 0.7:
                dialog.ListCalculated.item(i).setBackground(red07)
            elif 0.5 <=consv < 0.6:
                dialog.ListCalculated.item(i).setBackground(red06)
            elif 0.4 <=consv < 0.5:
                dialog.ListCalculated.item(i).setBackground(red05)
            elif 0.3 <=consv < 0.4:
                dialog.ListCalculated.item(i).setBackground(red04)
            elif 0.2 <=consv < 0.3:
                dialog.ListCalculated.item(i).setBackground(red03)
            elif 0.1 <=consv < 0.2:
                dialog.ListCalculated.item(i).setBackground(red02)
            else:
                dialog.ListCalculated.item(i).setBackground(red01)

    @staticmethod
    def recalculate_hetatm_clusters(worker):
        # Calculates all the clusters of heteratoms, run only when originaly running the superposition and when the epsilon value is changed
        global dict_hetatm_typ_coords_master
        global hetatm_typ_max_consv
        global hetatm_typ_nr
        global DBSCAN_eps_of_last_calc
        hetatm_typ_max_consv = {}
        hetatm_typ_nr = {}
        DBSCAN_eps_of_last_calc = dialog.SpinDB.value()
        report_list_1.append("\n-"*12 + "HETATM CLUSTERING"+"-"*12)
        report_list_1.append("Clustering of heteroatoms with DBSCAN, epslilon = {} A\n\n".format(DBSCAN_eps_of_last_calc))
        prog_count = 0
        worker.prog_sig.emit(0)
        worker.prog_msg.emit("Running clustering:")

        for hetatm_typ in dict_hetatm_typ_coords_master.keys():
            #print("starting",hetatm_typ)
            if hetatm_typ == "HOH-O":
                worker.prog_msg.emit("Running clustering for " + hetatm_typ+". Water may take a little longer")
            else:
                worker.prog_msg.emit("Running clustering for " + hetatm_typ)
                
            hetatm_typ_nonam, atom_name = hetatm_typ.split("-")

            array_coords = np.array(dict_hetatm_typ_coords_master[hetatm_typ])
            array_coords_x = array_coords[:,0]
            max_pop, max_consv =  HETATM_clustering.make_hetatm_iden_clusters(array_coords, hetatm_typ)
            
            
            
            if max_consv > 1:
                max_consv = 1


            if hetatm_typ_nonam not in hetatm_typ_max_consv:
                hetatm_typ_max_consv[hetatm_typ_nonam] = max_consv
                hetatm_typ_nr[hetatm_typ_nonam] = len(array_coords_x)
            else:
                hetatm_typ_nr[hetatm_typ_nonam] += len(array_coords_x)
                if hetatm_typ_max_consv[hetatm_typ_nonam] < max_consv:
                    hetatm_typ_max_consv[hetatm_typ_nonam] = max_consv
            prog_count += 1
            prog = int(prog_count/len(dict_hetatm_typ_coords_master.keys())*100)
            worker.prog_sig.emit(prog)
            if worker.exit_flag:
                
                break

        
        worker.finished.emit()
    
    @staticmethod
    def display_list_hetatm_types(dict_hetatm_typ_coords,display_calc_cl = True):
        # display the list of heteroatom types
        dialog.ListHetatms.clear()
        total_hetatms = 0

        for hetatm_typ in dict_hetatm_typ_coords.keys():
            total_hetatms += len(dict_hetatm_typ_coords[hetatm_typ])

        ion_list = Find.find_ion_list()
        nr_ions = 0


        hetams_clusters = []

        global hetatm_typ_max_consv
        global hetatm_typ_nr
            
        
        if display_calc_cl:
            for hetatm_typ in dict_hetatm_typ_coords.keys():
                HETATM_clustering.display_cluster_info_hetatm_typ(hetatm_typ)

        max_max_consv = 0
        max_ion_consv = 0

        for hetatm_typ, max_consv in hetatm_typ_max_consv.items():
            if hetatm_typ not in [x.split("-")[0] for x in dict_hetatm_typ_coords.keys()]:
                continue
            
            
            
            
            if max_consv == 0:
                red_nr = red_list[0]
            else:
                try:
                    red_nr = red_list[int(10*(max_consv+0.09999999))]
                except:
                    red_nr = red_list[0]

            text = '{}, {} present, max. consv.: {}'.format(hetatm_typ,hetatm_typ_nr[hetatm_typ],max_consv)
            
            hetams_clusters.append({'data':[text,red_nr],'max_consv':max_consv})
            if max_max_consv < max_consv:
                max_max_consv = max_consv
            if hetatm_typ.upper() in ion_list:
                nr_ions += hetatm_typ_nr[hetatm_typ]
                if max_ion_consv < max_consv:
                    max_ion_consv = max_consv


        

        dialog.ListHetatms.addItem('ALL HETATM, {} total, max. consv.: {}'.format(total_hetatms,max_max_consv))
        dialog.ListHetatms.addItem('IONS, {} present, max. consv.: {}'.format(nr_ions,max_ion_consv))

        dialog.ListHetatms.item(0).setBackground(red_list[int(10*(max_max_consv+0.09999999))])
        dialog.ListHetatms.item(1).setBackground(red_list[int(10*(max_ion_consv+0.09999999))])

        
        def return_max_consv(hetams_clusters):
            return hetams_clusters['max_consv']


        hetams_clusters.sort(reverse = True, key = return_max_consv)
        i = 2
        for dic in hetams_clusters:
            dialog.ListHetatms.addItem(dic['data'][0])
            dialog.ListHetatms.item(i).setBackground(dic['data'][1])
            i=i+1
    
    @staticmethod
    def run_Filter_HETATM():
        
        global popup_bar
        popup_bar.to_run = "Filter_HETATM"
        if DBSCAN_eps_of_last_calc != dialog.SpinDB.value():
            popup_bar.show_progress()

        popup_bar.start_progress()
        
    @staticmethod
    def Filter_HETATM_GUI_interaction_types():
        dict_hetatm_typ_coords = HETATM_clustering.reduce_dict_hetatm_typ_coords_master(True)
        HETATM_clustering.display_list_hetatm_types(dict_hetatm_typ_coords,False) # Re-display the heteroatom types, applying any changes (allowed/ not allowed list/ recalculation)
    
    
    @staticmethod
    def Filter_HETATM_GUI_interaction():
        try:
            selected_text = dialog.ListHetatms.currentItem().text()
        except AttributeError:
            selected_text = 'ALL HETATM, {} total, max. consv.: {}'
        selected_hetatm_typ = selected_text.split()[0][:-1]
        dict_hetatm_typ_coords = HETATM_clustering.reduce_dict_hetatm_typ_coords_master()

        dialog.ListCalculated.clear()
        #print("Filtering by: ",selected_text)
        
        if selected_text.split()[0] == 'ALL': # ALL 
            for hetatm_typ in dict_hetatm_typ_coords.keys():
                HETATM_clustering.display_cluster_info_hetatm_typ(hetatm_typ)
        elif selected_text.split()[0] == 'IONS,': # IONS
            ion_list = Find.find_ion_list()
            for hetatm_typ in dict_hetatm_typ_coords.keys():
                if hetatm_typ.split("-")[0].upper() in ion_list:
                    HETATM_clustering.display_cluster_info_hetatm_typ(hetatm_typ)
        else:

            for hetatm_typ in dict_hetatm_typ_coords.keys(): # specific heteroatom type
                if hetatm_typ.split("-")[0] == selected_hetatm_typ:
                    HETATM_clustering.display_cluster_info_hetatm_typ(hetatm_typ)



        #HETATM_clustering.display_list_hetatm_types(dict_hetatm_typ_coords,False) # Re-display the heteroatom types, applying any changes (allowed/ not allowed list/ recalculation)
        #if DBSCAN_eps_of_last_calc == dialog.SpinDB.value():
    @staticmethod
    def Filter_HETATM(worker):
        # Filter Heteroatoms button
        worker.prog_sig.emit(0)


        if DBSCAN_eps_of_last_calc != dialog.SpinDB.value(): # Recalculating heteroatom clusters, epsilon has changed
            print("Value of Epsilon parameter of DBSCAN changed\nRecalculating HETATM clusters")
            HETATM_clustering.recalculate_hetatm_clusters(worker)
        
        worker.Filter_HETATM_GUI.emit()
        worker.Filter_HETATM_GUI_types.emit()
        worker.finished.emit()

#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#
#                                                                                            Analysis
#
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################

class Analysis:
    @staticmethod
    def run_analyze_hetams():
        
        global popup_bar
        popup_bar.to_run = "analyze_hetams"
        popup_bar.show_progress()
        popup_bar.start_progress()


    @staticmethod
    def analyze_hetams(worker):
        # GO button, runs the superposition, clusters heteroatoms
        
        # cleanup of old files
        if platform == "win32":
            Run_Subprocess("DEL *.ent.gz /S")
            Run_Subprocess("DEL *.srf /S")
            Run_Subprocess("DEL *.chain.pdb /S")
            Run_Subprocess("DEL *.cons.pdb /S")
            Run_Subprocess("DEL Align_TM* /S")
            Run_Subprocess("DEL *.rota.pdb /S")
            Run_Subprocess("DEL AAA_NOSQL* /S")
            Run_Subprocess("DEL query.json /S")
            Run_Subprocess("DEL info.json /S")
            Run_Subprocess("DEL srfs.txt /S")
            Run_Subprocess("DEL *.stdout /S")
            Run_Subprocess("DEL Align_Deep* /S")

        else:
            Run_Subprocess("rm ./*.ent.gz")
            Run_Subprocess("rm ./*.srf")
            Run_Subprocess("rm ./*.chain.pdb")
            Run_Subprocess("rm ./*.cons.pdb")
            Run_Subprocess("rm ./Align_TM*")
            Run_Subprocess("rm ./*.rota.pdb")
            Run_Subprocess("rm ./AAA_NOSQL*")
            Run_Subprocess("rm ./query.json")
            Run_Subprocess("rm ./info.json")
            Run_Subprocess("rm ./srfs.txt")
            Run_Subprocess("rm ./*.stdout")
            Run_Subprocess("rm ./Align_Deep*")

        global dict_hetatm_typ_coords_master 
        dict_hetatm_typ_coords_master = {}
        global SELECTED_SITE
        global SELECTED_SITE_CHAIN
        global entities
        global report_list_1
        superposition_tool = dialog.ComboAlignTool.currentText()
        # NUM CPU:
        try:
            processors_available_local = str(multiprocessing.cpu_count())
        except:
            processors_available_local = "1"
        # NUM_CPUS //

        
        report_list_1.append( "\n\n" + "-" * 25 )
        report_list_1.append(Plugin_Name+" REPORT file")
        report_list_1.append("-" * 25 + "\n\n")
        
        
        machine = open (dl_machine_sett, "r")
        check_path = machine.read()
        machine.close()
        MHL_dir = os.path.join(check_path,"MHL")
        if not os.path.isdir(MHL_dir):
            os.mkdir(MHL_dir)
        report_dir = os.path.join(check_path,"Reports")
        if not os.path.isdir(report_dir):
            os.mkdir(report_dir)

        
        if report_list_1 == []:
            report_filename = "report_" + dialog.LineProtein.text().lower()+"_"+ dialog.ComboAlignTool.currentText()+ ".txt"
            nova_datoteka = open(os.path.join(report_dir,report_filename), "w").close()

        # single alligned chain per PDB entry
        one_or_multiple = dialog.CheckAligned.isChecked()
        # get setting on superposition starting point analysis
        bsite_space_check = dialog.CheckAnalyze.isChecked()

        chain_sel = dialog.CheckCompare.isChecked()
        

        target_complex_2 = dialog.LineProtein.text().lower()


        try:
            bsite_selection_full = dialog.ListIdentify.currentItem()
            bsite_selection = bsite_selection_full.text()
            chain_selection = bsite_selection_full.text().split(".")[-1]
            whole_chain_compare_selection = bsite_selection_full.text().split(" ")[0]
        except:
            worker.signal_make_msg_box.emit(Plugin_Name+" Warning"+ "|Invalid selection \n\nPlease select b-site or chain!")
            worker.finished.emit()
            return None

        if bsite_space_check == True:
            for element in BindingSites.bsite_unique_centers:
                if element[0] == bsite_selection:
                    SELECTED_SITE = element
                    SELECTED_SITE_CHAIN = str(chain_selection).upper()

        if bsite_space_check == False:
            for element in BindingSites.bsite_unique_centers:
                if element[0] == bsite_selection:
                    SELECTED_SITE = element
                    SELECTED_SITE_CHAIN = str(chain_selection).upper()
        
        if chain_sel == True:
            SELECTED_SITE = []
            SELECTED_SITE.append("no binding site used in analysis")
            SELECTED_SITE_CHAIN = str(whole_chain_compare_selection).upper()
            
        else:
            pass
        

        report_list_1.append("\n\n\nExamined complex: " + target_complex_2)
        report_list_1.append("Whole chain setting used: " + str(chain_sel))
        if chain_sel == True:
            report_list_1.append("Whole chain selection: " + whole_chain_compare_selection)
            report_list_1.append("Binding site selection: / (not used)")
            report_list_1.append("Chain selection: " + whole_chain_compare_selection)
        else:
            report_list_1.append("Whole chain selection: / (not used)")
            report_list_1.append("Binding site selection: " + bsite_selection)
            report_list_1.append("Chain selection: " + chain_selection)

        if  dialog.ComboMMseqs2.currentText() == "Custom Cluster":
            report_list_1.append("Used a custom clusters of complexes")
        else:
            report_list_1.append("Used PDB complexes with " + dialog.ComboMMseqs2.currentText() + " sequence identity")
        try:
            report_list_1.append("Unique structures in identified cluster: " + str(list(analyzesed_complexes_dict.keys())) + "\n\n\n")
        except NameError:
            worker.signal_make_msg_box.emit(Plugin_Name+" Warning"+ "|Click Find to define the cluster of complexes to be analysed")
            worker.finished.emit()
            return 
        
        if target_complex_2.lower() not in analyzesed_complexes_dict.keys():
            worker.signal_make_msg_box.emit(Plugin_Name+" Warning" "|Target complex not in defined cluster of complexes\nClick Find to redefine the cluster of complexes to be analysed")
            worker.finished.emit()
            return 

        for komp in analyzesed_complexes_dict.keys():
            try:
                #finding the identites of heteroatoms
                Find.find_hetnames(komp.lower() + ".pdb")
            except:
                pass


        global atom_max_x
        global atom_min_x
        global atom_max_y
        global atom_min_y
        global atom_max_z
        global atom_min_z

        bsite_rad = dialog.SpinAddBsiteRad.value()

        if chain_sel == True:
            correction_x = []
            correction_y = []
            correction_z = []
            for element in list_atoms_xyzchain:
                if SELECTED_SITE_CHAIN == element[3]:
                    correction_x.append(float(element[0]))
                    correction_y.append(float(element[1]))
                    correction_z.append(float(element[2]))
            atom_max_x = max(correction_x) + bsite_rad
            atom_min_x = min(correction_x) - bsite_rad
            atom_max_y = max(correction_y) + bsite_rad
            atom_min_y = min(correction_y) - bsite_rad
            atom_max_z = max(correction_z) + bsite_rad
            atom_min_z = min(correction_z) - bsite_rad
        else:
            atom_max_x = SELECTED_SITE[5] + bsite_rad
            atom_min_x = SELECTED_SITE[4] - bsite_rad
            atom_max_y = SELECTED_SITE[7] + bsite_rad
            atom_min_y = SELECTED_SITE[6] - bsite_rad
            atom_max_z = SELECTED_SITE[9] + bsite_rad
            atom_min_z = SELECTED_SITE[8] - bsite_rad





        
        use_curretn_MHL = dialog.CheckUseCurrentMHL.isChecked()

        timer_start() # Superpos
        #global popup_bar
        worker.prog_msg.emit("Running Superposition:")
        

        if not use_curretn_MHL:
            # running the superposition
            print(superposition_tool)
            if superposition_tool == "ProBiS":
                Align_methods.align_ProBiS(processors_available_local,bsite_selection,target_complex_2,chain_selection,analyzesed_complexes_dict,chain_sel,whole_chain_compare_selection,one_or_multiple,worker)
                print()
                print(Probis_cite)
                report_list_1.append("\n"+Probis_cite)
            elif superposition_tool[:5] == "PyMOL":
                print()

                #if dialog.CheckPyMolAlign.isChecked():
                if "align" in superposition_tool:
                    print(Pymol_align_cite)
                    report_list_1.append("\n"+Pymol_align_cite)
                else:
                    print(Pymol_super_cite)
                    report_list_1.append("\n"+Pymol_super_cite)
                Align_methods.align_Pymol(processors_available_local,bsite_selection,target_complex_2,chain_selection,analyzesed_complexes_dict,chain_sel,whole_chain_compare_selection,one_or_multiple,worker,superposition_tool)
            elif superposition_tool == "TM-align":
                print()
                print(TM_align_cite)
                report_list_1.append("\n"+TM_align_cite)
                Align_methods.align_TMalign(processors_available_local,bsite_selection,target_complex_2,chain_selection,analyzesed_complexes_dict,chain_sel,whole_chain_compare_selection,one_or_multiple,worker)
            elif superposition_tool == "GANGSTA+":
                print()
                print(gplus_cite)
                report_list_1.append("\n"+gplus_cite)
                Align_methods.align_gplus(processors_available_local,bsite_selection,target_complex_2,chain_selection,analyzesed_complexes_dict,chain_sel,whole_chain_compare_selection,one_or_multiple,worker)
                pml_filenames = glob("*.pml")
                del_str = ""
                for pml_filename in pml_filenames:
                    pdb_filename = pml_filename[:-3]+"pdb"
                    del_str += " {} {}".format(pml_filename,pdb_filename)
                if platform == "win32":
                    
                    Run_Subprocess("DEL"+del_str)

                else:
                    Run_Subprocess("rm"+del_str)
            elif superposition_tool == "DeepAlign":
                print()
                print(Deep_align_cite)
                report_list_1.append("\n"+Deep_align_cite)
                Align_methods.align_DeepAlign(processors_available_local,bsite_selection,target_complex_2,chain_selection,analyzesed_complexes_dict,chain_sel,whole_chain_compare_selection,one_or_multiple,worker)



            #######################################################################################
            open(os.path.join(MHL_dir,dialog.LineMHLName.text()),'w').close()
            def find_all_hetatms(filename): #for the .pdb file, that has B factors correct
                with open(os.path.join(MHL_dir,dialog.LineMHLName.text()),'a') as mhl_file:
                    mhl_file.write("#HETATMS in file {}:\n".format(filename))

                with open(filename,'r') as fil:
                    lines = fil.readlines()

                    vzorec_end_model = re.compile("^ENDMDL")

                    vzorec_hetatm =re.compile("^HETATM",re.IGNORECASE) 
                    for lin in lines:
                        match_end = vzorec_end_model.match(lin)
                        if match_end: 
                            break

                        match_hetatm = vzorec_hetatm.match(lin)
                        if match_hetatm:
                            typ = lin[17:20].strip()
                            x_hetm = float(lin[30:38])
                            y_hetm = float(lin[38:46])
                            z_hetm = float(lin[46:54])
                            B_hetm = float(lin[60:66])
                            seq_nr = int(lin[22:26])
                            ser_nr = int(lin[6:11])
                            chain_iden = lin[21:22].strip()
                            atom_name  = lin[12:16].strip()
                            with open(os.path.join(MHL_dir,dialog.LineMHLName.text()),'a') as mhl_file:
                                mhl_file.write('{} {} {} {} {} {} {} {} {} {}\n'.format(x_hetm,y_hetm,z_hetm,B_hetm,ser_nr,chain_iden,typ,atom_name,seq_nr,filename[:-4]))
                            HETATM_clustering.add_to_dict_hetatm_typ_coords_master(x_hetm,y_hetm,z_hetm,B_hetm,ser_nr,chain_iden,typ,atom_name,seq_nr,filename[:-4])
                            


            def find_all_hetatms_B_corr(filename_rota_pdb,filename_pdb,vzorec_end): # for the .rota.pdb files that have B factors 0, it takes B factors from filename_pdb
                with open(os.path.join(MHL_dir,dialog.LineMHLName.text()),'a') as mhl_file:
                    mhl_file.write("#HETATMS in file {}:\n".format(filename_rota_pdb))
                fil_rota =  open(filename_rota_pdb,'r')
                lines_rota = fil_rota.readlines()

                fil_pdb = open(filename_pdb,'r')
                lines_pdb_temp = fil_pdb.readlines()
                vzorec_hetatm =re.compile("^HETATM",re.IGNORECASE) 
                lines_pdb = []

                pdb_hetatm_dict = {}
                for lin_nr_pdb in range(len(lines_pdb_temp)): #search thorugh .pdb file
                    lin_pdb =  lines_pdb_temp[lin_nr_pdb]
                    match_hetatm_pdb = vzorec_hetatm.match(lin_pdb)
                    if match_hetatm_pdb:
                        lines_pdb.append(lin_pdb)
                        atom_name  = lin_pdb[12:16].strip()

                        typ = lin_pdb[17:20].strip() 
                        seq_nr = int(lin_pdb[22:26])
                        chain_iden = lin_pdb[21:22].strip()
                        pdb_hetatm_dict[(typ,seq_nr,atom_name)] = lin_pdb


                vzorec_end_model = re.compile(vzorec_end)

                fil_pdb.close()
                fil_rota.close()

                for lin_rota in lines_rota:
                    match_end = vzorec_end_model.match(lin_rota)
                    if match_end: 
                        break

                    match_hetatm = vzorec_hetatm.match(lin_rota)
                    if match_hetatm:                            #find hetatm in .rota.pdb file
                        typ = lin_rota[17:20].strip()
                        x_hetm = float(lin_rota[30:38])
                        y_hetm = float(lin_rota[38:46])
                        z_hetm = float(lin_rota[46:54])
                        B_hetm = 0
                        seq_nr = int(lin_rota[22:26])
                        ser_nr = int(lin_rota[6:11])
                        chain_iden = lin_rota[21:22].strip()
                        atom_name  = lin_rota[12:16].strip()


                        if (typ,seq_nr,atom_name) in pdb_hetatm_dict: #check if hetatms in .pdb and .rota.pdb are the same
                            
                            lin_pdb = pdb_hetatm_dict[(typ,seq_nr,atom_name)]
                            B_hetm = float(lin_pdb[60:66])
                            ser_nr = int(lin_pdb[6:11])
                            chain_iden = lin_pdb[21:22].strip()
                            
                        with open(os.path.join(MHL_dir,dialog.LineMHLName.text()),'a') as mhl_file:
                            #x y z B serial_nr chain_iden hetatom_typ atom_name sequence_nr PDB_id
                            mhl_file.write('{} {} {} {} {} {} {} {} {} {}\n'.format(x_hetm,y_hetm,z_hetm,B_hetm,ser_nr,chain_iden,typ,atom_name,seq_nr,filename_pdb[:-4]))
                            HETATM_clustering.add_to_dict_hetatm_typ_coords_master(x_hetm,y_hetm,z_hetm,B_hetm,ser_nr,chain_iden,typ,atom_name,seq_nr,filename_pdb[:-4])

                
                    

            
            filenames = glob("*.rota.pdb")
            # the last one is the target complex
            vzorec_end ="^ENDMDL"
            if superposition_tool == "TM-align":
                vzorec_end =  "^TER" # TM-aling has TER instead od ENDMDL at the end

            print("Alignment of complexes complete, finding HETATMS in aligned complexes:")
            i = 1
            prog_count = 0
            worker.prog_sig.emit(prog_count)
            worker.prog_msg.emit("Alignment of complexes complete, finding HETATMS in aligned complexes:")
            for filename_rota in filenames: 
                print("  Finding heteratoms in {}, {}/{}".format(filename_rota,i,int(len(set(filenames)))+1))
                # finding heteroatoms
                pdb_fnm = filename_rota.split("_")[-1].split(".")[0][:-1]+".pdb"
                find_all_hetatms_B_corr(filename_rota,pdb_fnm,vzorec_end) # takes B from the original pdb file, not the aligned one
                i += 1
                prog_count += 1
                prog = int(prog_count/(len(filenames)+1)*100)
                worker.prog_sig.emit(prog)
                if worker.exit_flag:
                    
                    break



            filenames.append(target_complex_2 + ".pdb") # the target complex
            print("  Finding heteratoms in {}, {}/{}".format(target_complex_2 + ".pdb",i,int(len(set(filenames)))))
            prog_count += 1
            prog = int(prog_count/(len(filenames)+1)*100)
            worker.prog_sig.emit(prog)
            if worker.exit_flag:
                
                return


            find_all_hetatms(target_complex_2 + ".pdb") # finds heteroatoms in the target complex
            entities = int(len(set(filenames)))
            print("Studied filenames: ",filenames)
            print("Number of collected structures: ",entities)


            ##########################################################################################

        else: # no roatiton, reading form Master Heteroatom List
            no_rot_cite = "No alignment of complexes\nReading data from Master Heteroatom List:  {}".format(dialog.LineMHLName.text())
            print()
            print(no_rot_cite)
            report_list_1.append("\n"+no_rot_cite)
            unique_files = []
            with open(os.path.join(MHL_dir,dialog.LineMHLName.text()),'r') as mhl_file:
                lines = mhl_file.readlines()
                
                for line in lines:
                    if line[0] == "#":                    #comments in MHL; #HETATMS in file: {}".format(filename_rota_pdb)
                        filename_rota_pdb = line.strip().split(" ")[3]
                        if filename_rota_pdb not in unique_files:
                            unique_files.append(filename_rota_pdb)
            entities = len(unique_files)
            dict_hetatm_typ_coords_master = Find.find_hetams_from_master_list()

        print("Protein superposition took {} s".format(time_since_timer_start()))
        timer_start()

        # bounding box object, displayed with box button
        global boundingBox
        boundingBox = [LINEWIDTH, 2.0, BEGIN, LINES,
                COLOR, float(1), float(0), float(0),

                VERTEX, atom_min_x, atom_min_y, atom_min_z,       #1
                VERTEX, atom_min_x, atom_min_y, atom_max_z,       #2

                VERTEX, atom_min_x, atom_max_y, atom_min_z,       #3
                VERTEX, atom_min_x, atom_max_y, atom_max_z,       #4

                VERTEX, atom_max_x, atom_min_y, atom_min_z,       #5
                VERTEX, atom_max_x, atom_min_y, atom_max_z,       #6

                VERTEX, atom_max_x, atom_max_y, atom_min_z,       #7
                VERTEX, atom_max_x, atom_max_y, atom_max_z,       #8


                VERTEX, atom_min_x, atom_min_y, atom_min_z,       #1
                VERTEX, atom_max_x, atom_min_y, atom_min_z,       #5

                VERTEX, atom_min_x, atom_max_y, atom_min_z,       #3
                VERTEX, atom_max_x, atom_max_y, atom_min_z,       #7

                VERTEX, atom_min_x, atom_max_y, atom_max_z,       #4
                VERTEX, atom_max_x, atom_max_y, atom_max_z,       #8

                VERTEX, atom_min_x, atom_min_y, atom_max_z,       #2
                VERTEX, atom_max_x, atom_min_y, atom_max_z,       #6


                VERTEX, atom_min_x, atom_min_y, atom_min_z,       #1
                VERTEX, atom_min_x, atom_max_y, atom_min_z,       #3

                VERTEX, atom_max_x, atom_min_y, atom_min_z,       #5
                VERTEX, atom_max_x, atom_max_y, atom_min_z,       #7

                VERTEX, atom_min_x, atom_min_y, atom_max_z,       #2
                VERTEX, atom_min_x, atom_max_y, atom_max_z,       #4

                VERTEX, atom_max_x, atom_min_y, atom_max_z,       #6
                VERTEX, atom_max_x, atom_max_y, atom_max_z,       #8

                END
        ]


        dx = abs(atom_max_x - atom_min_x)
        dy = abs(atom_max_y - atom_min_y)
        dz = abs(atom_max_z - atom_min_z)
        system_volume = dx * dy * dz
        #######################################################################################################################################
        #######################################################################################################################################
        #######################################################################################################################################
        #######################################################################################################################################
        #######################################################################################################################################
        #######################################################################################################################################

        
        report_list_1.append("System volume is: %d cubic A\n" % (system_volume))




        ########################################################################################################
        # ----------------------------------------------------------------------
        
        HETATM_clustering.recalculate_hetatm_clusters(worker)        

        worker.analysis_GUI.emit([target_complex_2,chain_sel,whole_chain_compare_selection,bsite_selection,chain_selection,superposition_tool,system_volume])
        
        return None

    def Analysis_GUI_interaction(analysis_GUI):

        args_ls = list(analysis_GUI)

        target_complex_2,chain_sel,whole_chain_compare_selection,bsite_selection,chain_selection,superposition_tool,system_volume = args_ls
        dialog.PlainInfo.clear() #######################################################################################################################################
        dialog.Tabs.setCurrentIndex(1) #######################################################################################################################################

        dialog.PlainInfo.insertPlainText("Examined complex: " + target_complex_2+"\n")
        dialog.PlainInfo.insertPlainText("Whole chain setting used: " + str(chain_sel)+"\n")
        if chain_sel == True:
            dialog.PlainInfo.insertPlainText("Whole chain selection: " + whole_chain_compare_selection+"\n")
            dialog.PlainInfo.insertPlainText("Binding site selection: / (not used)"+"\n")
            dialog.PlainInfo.insertPlainText("Chain selection: " + whole_chain_compare_selection+"\n")
        else:
            dialog.PlainInfo.insertPlainText("Whole chain selection: / (not used)"+"\n")
            dialog.PlainInfo.insertPlainText("Binding site selection: " + bsite_selection+"\n")
            dialog.PlainInfo.insertPlainText("Chain selection: " + chain_selection+"\n")

        dialog.PlainInfo.insertPlainText("Superposition method used: {}\n".format(superposition_tool))
        if  dialog.ComboMMseqs2.currentText() == "Custom Cluster":
            dialog.PlainInfo.insertPlainText("Used a custom clusters of complexes\n")
        else:
            dialog.PlainInfo.insertPlainText("Used PDB complexes with " + dialog.ComboMMseqs2.currentText() + " sequence identity\n")

        #dialog.PlainInfo.insertPlainText("Used PDB clusters with: " + dialog.ComboMMseqs2.currentText() + " %"+"\n")
        dialog.PlainInfo.insertPlainText("Unique structures in identified cluster: " + str(list(analyzesed_complexes_dict.keys())) + "\n\n")
        dialog.PlainInfo.insertPlainText("System volume is: %d cubic A\n" % (system_volume))

        dialog.ListCalculated.clear()
        

        HETATM_clustering.display_list_hetatm_types(HETATM_clustering.reduce_dict_hetatm_typ_coords_master(True),display_calc_cl = True)
        print("HETATM clustering took {} s".format(time_since_timer_start()))


        list_all_hetatms_types = list(set(map(lambda x: x.split("-")[0],dict_hetatm_typ_coords_master.keys())))
        dialog.PlainInfo.insertPlainText("\n\nAll HETATM types present: {}\n".format(list_all_hetatms_types))

    
        try:
            # print the identites of hetaroatoms
            dialog.PlainInfo.insertPlainText("\nIdentities of identified heteroatoms:\n")
            for typ in dic_all_hetatm_types_names.keys():
                dialog.PlainInfo.insertPlainText("HETATM Type:  {}\n   Name:   {}\n   Formula:    {}\n".format(typ,dic_all_hetatm_types_names[typ],dic_all_hetatm_types_formulas[typ]))
        except:
            pass
        dialog.Tabs.setCurrentIndex(1)

        ################################################################################
    
        
        
        

    





#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#
#                                                                                            PyMOL Interface
#
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################

colour_counter = 0


class pyMOLinterface:
    """use wonderful PyMOL for visualisation of collected results"""
    """thanks! Warren L. DeLano!"""

    @staticmethod
    def PyMOL_close_resi_contacts(): 
        # draws distances to the closest atoms of close residues
        target_complex_3 = dialog.LineProtein.text().lower()
        cmd.delete("dist*")

        def dist(coords1,coords2):
            return np.sqrt((coords1[0]-coords2[0])**2+(coords1[1]-coords2[1])**2+(coords1[2]-coords2[2])**2)

        stored.list_clus = []

        cmd.iterate_state("1","clus*","stored.list_clus.append([x,y,z,ID,index,model  ])")
        for x_cl,y_cl,z_cl,ID_cl,index_cl,model_cl in stored.list_clus:
            stored.list = []
            cmd.select("bsite", "((model {} and id {} and index {}) around 6) and (polymer or organic) and {}".format(model_cl,ID_cl,index_cl,target_complex_3))
            cmd.iterate_state("1","bsite","stored.list.append([x,y,z,resi,ID,index])")
            cmd.delete("bsite")


            resi_dic = {}
            for x,y,z,resi,ID,index in stored.list:
                if resi in resi_dic:

                    resi_dic[resi].append([x,y,z,resi,ID,index])
                else:

                    resi_dic[resi] = [[x,y,z,resi,ID,index]]
            min_dist_total_list = []

            for resi in resi_dic.keys():
                min_dist = 10000
                min_list = []
                for x,y,z,resi,ID,index in resi_dic[resi]:
                    ds = dist([x_cl,y_cl,z_cl],[x,y,z])
                    if ds < min_dist:
                        min_dist = ds
                        min_list = [x,y,z,resi,ID,index]
                
                min_dist_total_list.append(min_list)
            
            


            sele_cmd = ""
            for i in range(len(min_dist_total_list)):
                resi  = min_dist_total_list[i][-3]
                ID =  min_dist_total_list[i][-2]
                index = min_dist_total_list[i][-1]
                sele_cmd +="(id {} and index {}) ".format(ID,index)
                if i != len(min_dist_total_list)-1:
                    sele_cmd += "or "
            
            cmd.select("bsites_closest_atom", sele_cmd)
            cmd.deselect()

            cmd.distance("dist_closest_bsite_atoms","bsites_closest_atom","model {} and id {} and index {}".format(model_cl,ID_cl,index_cl),"4")

        cmd.distance("dist_inter_cluts", "clus*","clus*", "4.0")

        cmd.set("dash_color", "magenta")
        cmd.set("dash_gap", "0.2")
        cmd.set("dash_length", "0.2")
        cmd.set("dash_round_ends","on")
        cmd.set("dash_width","3")


    @staticmethod
    def pyMOL_fetch_system():
        # deletes all, fetches the targed complex

        cmd.delete(name = "all")
        target_complex_3 = dialog.LineProtein.text().lower()
        print("Target Complex: "+target_complex_3)

        if custom == True:
            cmd.load(filename, target_complex_3)
        else:
            cmd.fetch(target_complex_3, target_complex_3)
        
        cmd.show ("cartoon", target_complex_3)        
        cmd.set("cartoon_color", "white", target_complex_3)
        cmd.hide("lines", "all")
        cmd.util.cbag(selection = target_complex_3)
        cmd.show("surface", target_complex_3)
        cmd.set("transparency", "0.9")
        cmd.set("surface_color", "white")
        cmd.show("sticks", "organic")
        cmd.color("blue", "organic")
        waters = "{}_waters".format(target_complex_3)
        cmd.select(waters, "resn hoh")
        cmd.show("nonbonded", waters)
        cmd.deselect()


    @staticmethod
    def pyMOL_chain_box():
        # loads the bounding box
        cmd.load_cgo(boundingBox, "box")

    @staticmethod
    def pyMOL_bsite_cluster():
        # display binding sites of clusters
        cmd.select("bsites", "clus* around 6")
        cmd.select("byres bsites")
        cmd.show("sticks", "byres bsites")
        cmd.util.cbay("byres bsites")
        cmd.set_bond("stick_radius", "0.1", "byres bsites")
        cmd.select("sele", "name ca and byres bsites")
        cmd.show("sticks", "organic")
        cmd.color("blue", "organic")
        cmd.util.cnc ("organic")
        cmd.set_bond("stick_radius", "0.25", "organic")
        
    @staticmethod
    def pyMOL_display_cluster():
        # displays cluster(s) and writes report
        display_clusters_setting = dialog.CheckKeep.isChecked()
        bsite_space_check = dialog.CheckAnalyze.isChecked()
        chain_sel = dialog.CheckCompare.isChecked()

        # B factor
        debye_waller_check = dialog.CheckDebye.isChecked()

        master_bsite_list_het = []
        master_bsite_list_het_coord_x = []
        master_bsite_list_het_coord_y = []
        master_bsite_list_het_coord_z = []
        master_list_het = []
        master_list_het_coord_x = []
        master_list_het_coord_y = []
        master_list_het_coord_z = []
        # Debye Waller
        master_bsite_list_atom_iso_displacement = []
        master_list_het_iso_disp = []
        # master_list_names = []
        master_bsite_list_info = []
        master_list_info = []


        



        try:
            selected_hetatm_typ = dialog.ListCalculated.currentItem().text().split()[4]
           
            if selected_hetatm_typ == "least":  selected_hetatm_typ = dialog.ListCalculated.currentItem().text().split()[6]

        except AttributeError:
            QtWidgets.QMessageBox.about(dialog, Plugin_Name + " WARNING", "Select HETATM cluster(s) to display in the PyMOL viewer")
            return
        print('Clusters of {} Heteroatoms:'.format(selected_hetatm_typ))

        mhl_file = open(os.path.join(MHL_dir,dialog.LineMHLName.text()), "r") 
        lines = mhl_file.readlines()
                #f'{x_hetm} {y_hetm} {z_hetm} {B_hetm} {filename} {seq_nr} {typ}\n'
        
        vdw_rad = dialog.SpinRadius.value()
        for line in lines:
            if line[0] == "#":
                continue
            x_hetm,y_hetm,z_hetm,B_hetm,serial_nr,chain_iden,hetatm_typ,atom_name,seq_nr,PDB_id = line.strip().split(" ")
                            #x y z B serial_nr chain_iden hetatom_typ atom_name sequence_nr PDB_id
            temp_list = []
            if hetatm_typ+"-"+atom_name == selected_hetatm_typ:
            
                x = float(x_hetm)
                y = float(y_hetm)
                z = float(z_hetm)
                B = float(B_hetm)
                
                if B < 0:
                    B = 0
                info = "{} {} {} {} {} {}".format(serial_nr,chain_iden,hetatm_typ,atom_name,seq_nr,PDB_id)


                
                isotropni_displacement = math.sqrt(B/(8*((math.pi)**2))) + vdw_rad
                # B = (isotropni_displacement-vdw_rad)**2 * (8*((math.pi)**2)



                temp_list.append(x)
                temp_list.append(y)
                temp_list.append(z)

                # name of bsite, axerage x, average y, average z, min x, max x, min y, max y, min z, max z


                if bsite_space_check == True and chain_sel == False:
                    #if SELECTED_SITE[4] - 4 <= temp_list[0] <= SELECTED_SITE[5] + 4:
                    #    if SELECTED_SITE[6] - 4 <= temp_list[1] <= SELECTED_SITE[7] + 4:
                    #        if SELECTED_SITE[8] - 4 <= temp_list[2] <= SELECTED_SITE[9] + 4:
                    if atom_min_x <= temp_list[0] <= atom_max_x:
                        if atom_min_y <= temp_list[1] <= atom_max_y:
                            if atom_min_z <= temp_list[2] <= atom_max_z:
                                master_bsite_list_het.append(temp_list)
                                x2 = float(temp_list[0])
                                master_bsite_list_het_coord_x.append(x2)
                                y2 = float(temp_list[1])
                                master_bsite_list_het_coord_y.append(y2)
                                z2 = float(temp_list[2])
                                master_bsite_list_het_coord_z.append(z2)
                                master_bsite_list_atom_iso_displacement.append(isotropni_displacement)
                                master_bsite_list_info.append(info)

                if bsite_space_check == False and chain_sel == False:
                    if atom_min_x <= temp_list[0] <= atom_max_x:
                        if atom_min_y <= temp_list[1] <= atom_max_y:
                            if atom_min_z <= temp_list[2] <= atom_max_z:
                                master_list_het_coord_x.append(x)
                                master_list_het_coord_y.append(y)
                                master_list_het_coord_z.append(z)
                                master_list_het.append(temp_list)
                                master_list_het_iso_disp.append(isotropni_displacement)
                                master_list_info.append(info)

                if chain_sel == True:
                    if atom_min_x <= temp_list[0] <= atom_max_x:
                        if atom_min_y <= temp_list[1] <= atom_max_y:
                            if atom_min_z <= temp_list[2] <= atom_max_z:
                                master_list_het_coord_x.append(x)
                                master_list_het_coord_y.append(y)
                                master_list_het_coord_z.append(z)
                                master_list_het.append(temp_list)
                                master_list_het_iso_disp.append(isotropni_displacement)
                                master_list_info.append(info)


                else:
                    pass


        if bsite_space_check == True and chain_sel == False:
            master_list_het = master_bsite_list_het
            master_list_het_iso_disp = master_bsite_list_atom_iso_displacement
            master_list_info = master_bsite_list_info
        else:
            pass

        mhl_file.close()


        try:
            # example:
            # 16 clusters with 16 HOH-O hetatms. consv. 0.94
            try:
                cluster_selection = int(dialog.ListCalculated.currentItem().text().split()[3])
            except:
                cluster_selection = int(dialog.ListCalculated.currentItem().text().split()[5]) # at least
            try:
                consv_of_cluster = float(dialog.ListCalculated.currentItem().text().split()[7])
            except:
                consv_of_cluster = float(dialog.ListCalculated.currentItem().text().split()[9]) # at least

            print(cluster_selection,consv_of_cluster)
            # report list
            bsite_rad = dialog.SpinAddBsiteRad.value()
            report_list_1.append("\nBinding site info (name, avg x, y, z, min x, max x, min y, max y, min z, max z; box {} A around extremes): \n".format(bsite_rad) + str(SELECTED_SITE))
            report_list_1.append("\nExamined clusters with " + str(cluster_selection) + " or more {} heteroatoms\n".format(selected_hetatm_typ))
            report_list_1.append("-" * 25)

        except:
            QtWidgets.QMessageBox.about(dialog, Plugin_Name+" Warning", "Please select clusters to display")
            return
        selected_eps=dialog.SpinDB.value()
        print(selected_eps,cluster_selection)
        print("master_list_het",master_list_het)
        labels3D = DBSCAN(eps=selected_eps, min_samples=cluster_selection).fit_predict(np.array(master_list_het))

        i = 0
        tocke = []
        for element in labels3D:
            temp = []
            if element != -1:
                temp.append(master_list_het[i])
                temp.append(element)
                temp.append(master_list_het_iso_disp[i])
                temp.append(master_list_info[i])
                # tocke -> [x,y,z], nr_cl, B, info
                tocke.append(temp)

            else:
                pass
            i += 1


            
        global colour_counter
        
        if display_clusters_setting == False:
            colour_counter = 0
        
            
        
        def get_colour(consv):
            global colour_counter
            
            col_ls = ["0xFF||","0xFFFF|","0xFF|FF","0x|FF|","0x|FFFF","0x||FF"]
            
            colour_hex = hex(int(((1-float(consv))*255))).replace("0x","").upper()
            if len(colour_hex) == 1:
                colour_hex = "0" + colour_hex
            colour = col_ls[colour_counter%len(col_ls)].replace("|",colour_hex)
            
            
            
            return colour


        if debye_waller_check == False:

            if display_clusters_setting == False:
                cmd.delete("clus*")
                cmd.delete("iso_disp")

            else:
                pass

            for element in list(set(labels3D)):
                cluster_temp = []
                
                
                for sub_element in tocke:
                    if sub_element[1] == element:
                        cluster_temp.append(sub_element[0])
                        


                cluster_temp = np.array(cluster_temp)
                
                if int(element) == -1:
                    continue

                try:

                    print("Cluster {}, nr. of hetatm in cluster: {}".format(element,len(cluster_temp)))
                    print(cluster_temp)
                    clust_avg = [np.average(cluster_temp[:,0]),np.average(cluster_temp[:,1]),np.average(cluster_temp[:,2])]
                    clust_st_dev = [np.std(cluster_temp[:,0]),np.std(cluster_temp[:,1]),np.std(cluster_temp[:,2])]
                    print("Cluster {} averege pos:   {} {} {}".format(element,round(clust_avg[0],3),round(clust_avg[1],3),round(clust_avg[2],3)))
                    print("Cluster {} averege stdev: {} {} {}".format(element,round(clust_st_dev[0],3),round(clust_st_dev[1],3),round(clust_st_dev[2],3)))
                    if element != -1:
                        report_list_1.append("\n\nCluster nr. {}, composed of {} {} heteroatoms, conservation {}:".format(element,len(cluster_temp),selected_hetatm_typ,round(float(len(cluster_temp))/float(entities),3 ) ))
                        report_list_1.append("(x y z iso_disp serial_number chain_identifier hetatm_type hetatom_name sequence_number PDB_id)\n")
                        for sub_element in tocke:
                            if sub_element[1] == element:
                                # tocke -> [x,y,z], nr_cl, B, info
                                report_line = "{} {} {} {} {}".format(sub_element[0][0],sub_element[0][1],sub_element[0][2],round((float(sub_element[2])-vdw_rad)**2 * 8*(math.pi)**2,2),sub_element[3])
                                report_list_1.append(report_line)
                    
                    report_list_1.append("\nCluster {}   Averege pos.:   {} {} {}".format(element,round(clust_avg[0],3),round(clust_avg[1],3),round(clust_avg[2],3)))
                    report_list_1.append("Cluster {} Averege st. dev.: {} {} {}".format(element,round(clust_st_dev[0],3),round(clust_st_dev[1],3),round(clust_st_dev[2],3)))

                    #cmd.set_color("clus_color", "[%f, %f, %f]" % get_colour(round(float(len(cluster_temp))/float(entities))))#(1.0, (1.0 - consv_of_cluster), (1.0 - consv_of_cluster)))
                    colour = get_colour(round(float(len(cluster_temp))/float(entities),3))
                    pymol.cmd.do("pseudoatom clus_%s-%d_%.2f, vdw=%f, color=%s, pos=[%f, %f, %f]" % (selected_hetatm_typ, element, consv_of_cluster,vdw_rad,colour, clust_avg[0],clust_avg[1],clust_avg[2]))
                    cmd.show("spheres", "clus_%s-%d*" % (selected_hetatm_typ,element))

                except IndexError:
                    pass

        else:
            if display_clusters_setting == False:
                cmd.delete("clus*")
                cmd.delete("iso_disp")
            
            for element in list(set(labels3D)):
                cluster_temp = []
                for sub_element in tocke:
                    if sub_element[1] == element:
                        sub_element[0].append(sub_element[2])
                        cluster_temp.append(sub_element[0])

                if int(element) == -1:
                    continue

                cluster_temp = np.array(cluster_temp)
                try:
                    print("Cluster {}, nr. of hetatm in cluster: {}".format(element,len(cluster_temp)))
                    print(cluster_temp)
                    clust_avg = [np.average(cluster_temp[:,0]),np.average(cluster_temp[:,1]),np.average(cluster_temp[:,2])]
                    clust_st_dev = [np.std(cluster_temp[:,0]),np.std(cluster_temp[:,1]),np.std(cluster_temp[:,2])]


                    if element != -1:
                        report_list_1.append("\n\nCluster nr. {}, composed of {} {} heteroatoms, conservation {}:".format(element,len(cluster_temp),selected_hetatm_typ,round(float(len(cluster_temp))/float(entities),3 ) ))
                        report_list_1.append("(x y z iso_disp serial_number chain_identifier hetatm_type hetatom_name sequence_number PDB_id)\n")
                        for sub_element in tocke:
                            if sub_element[1] == element:
                                # tocke -> [x,y,z], nr_cl, B, info

                                report_line = "{} {} {} {} {}".format(sub_element[0][0],sub_element[0][1],sub_element[0][2],round((float(sub_element[2])-vdw_rad)**2 * 8*(math.pi)**2,2),sub_element[3])
                                report_list_1.append(report_line)
                    
                    print("Cluster {} averege pos:   {} {} {}".format(element,round(clust_avg[0],3),round(clust_avg[1],3),round(clust_avg[2],3)))
                    print("Cluster {} averege stdev: {} {} {}".format(element,round(clust_st_dev[0],3),round(clust_st_dev[1],3),round(clust_st_dev[2],3)))
                    report_list_1.append("\nCluster {}   Averege pos.:   {} {} {}".format(element,round(clust_avg[0],3),round(clust_avg[1],3),round(clust_avg[2],3)))
                    report_list_1.append("Cluster {} Averege st. dev.: {} {} {}".format(element,round(clust_st_dev[0],3),round(clust_st_dev[1],3),round(clust_st_dev[2],3)))
                    #cmd.set_color("clus_color", "[%f, %f, %f]" % (1.0, (1.0 - consv_of_cluster), (1.0 - consv_of_cluster)))
                    colour = get_colour(round(float(len(cluster_temp))/float(entities),3))
                    
                    pymol.cmd.do("pseudoatom clus_%s-%d_%.2f, vdw=%f, color=%s, pos=[%f, %f, %f]" % (selected_hetatm_typ, element, consv_of_cluster,vdw_rad,colour, clust_avg[0],clust_avg[1],clust_avg[2]))
                    cmd.show("spheres", "clus_%s-%d*" % (selected_hetatm_typ,element))

                    for tocka in cluster_temp:

                        pymol.cmd.do("pseudoatom iso_disp, vdw=%f, color=red, pos=[%f, %f, %f]" % (tocka[3], tocka[0],tocka[1],tocka[2]))

                    cmd.show("dots", "iso_disp")

                except IndexError:
                    pass
        colour_counter += 1
        

        # report on cluster
        report_filename = "report_" + dialog.LineProtein.text().lower()+"_"+ dialog.ComboAlignTool.currentText()+ ".txt"
        nova_datoteka = open(os.path.join(report_dir,report_filename), "w")
        for linija in report_list_1:
            nova_datoteka.write("%s\n" % linija)
        nova_datoteka.close()
        print("report created...")

