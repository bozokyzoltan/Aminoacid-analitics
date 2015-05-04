#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2015.04.17.

"""

from aa_data import AminoacidAnalyzes_file
from aa_data import Patient
import numpy as np
from excel import Excel
import wx
import os
import glob


class FrontEnd(wx.Frame):
    """
    """
    def __init__(self):
        """
        """

        wx.Frame.__init__(self, None, wx.ID_ANY, "Aminoacid analyzes",
                          size = (200, 80), style = wx.CAPTION | wx.MINIMIZE_BOX |wx.CLOSE_BOX)
        self.Centre()
        #
        self.panel1 = wx.Panel(self, wx.ID_ANY, size =(300,300))
        # ---------------------
        # 1. Load datafile
        # ---------------------
        self.button_data = wx.Button(self.panel1, id=wx.ID_ANY, label="Select data folder", pos=(5,5), size=(150,30))
        self.button_data.Bind(wx.EVT_BUTTON, self.onButton_data)
        #
        wx.StaticText(self.panel1, id=wx.ID_ANY, pos=(155,60), label='bozoky')
        #
        return None
    ### ========================================================================
    def onButton_data(self, event):
        """
        This method is fired when its corresponding button is pressed
        """
        openDirDialog = wx.DirDialog(self.panel1, 'Choose input directory', '.', wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        #
        if openDirDialog.ShowModal() != wx.ID_CANCEL:
            # ----------------------------------------
            # Get all data files
            # ----------------------------------------
            patients = {}
            for filename in glob.glob(openDirDialog.GetPath() + os.sep + '*.txt'):
                file_content = AminoacidAnalyzes_file(filename)
                if file_content.patient not in patients:
                    patients[file_content.patient] = Patient()
                patients[file_content.patient].add_data(file_content)
            #
            # ----------------------------------------
            # Create result folder
            # ----------------------------------------
            folder = ''.join((openDirDialog.GetPath(), os.sep, '..', os.sep, 'result', os.sep))
            if not os.path.isdir(folder):
                os.makedirs(folder)
            # ----------------------------------------
            # Go through each patient data and create an excel file for each
            # ----------------------------------------
            uptakes = {}
            for patient in patients:
                print patient
                patients[patient].calculate_average()
                # ----------------------------------------
                # Create Excel file
                # ----------------------------------------
                excelfile = Excel()
                # ----------------------------------------
                # Raw data
                # ----------------------------------------
                excelfile.add_sheet('raw data', patients[patient].raw_data()[1])
                # ----------------------------------------
                # Concentration
                # ----------------------------------------
                excelfile.add_sheet('concentration', patients[patient].raw_data()[0])
                # ----------------------------------------
                # Mean + Std
                # ----------------------------------------
                excelfile.add_sheet('mean', patients[patient].mean_std())
                # ----------------------------------------
                # Normalized to aminoacid
                # ----------------------------------------
                for i in xrange(len(patients[patient].aa)):
                    #
                    new = np.array(patients[patient].normalized(-1, i))
                    if i == 0:
                        a = new
                    else:
                        a = np.concatenate((a, new), axis = 1)
                #
                excelfile.add_sheet('normalized', a)
                # ----------------------------------------
                # Uptake
                # ----------------------------------------
                excelfile.add_sheet('uptake', patients[patient].normalized(0, -1))
                # ----------------------------------------
                # Uptake normalized
                # ----------------------------------------
                for i in xrange(len(patients[patient].aa)):
                    new = np.array(patients[patient].normalized(0, i))
                    if i == 0:
                        a = new
                    else:
                        a = np.concatenate((a, new), axis = 1)
                #
                excelfile.add_sheet('uptake normalized', a)
                # ----------------------------------------
                # Save Excel datafile
                # ----------------------------------------
                excelfile.save_excel_file(folder + patient + '.xls')
                # ----------------------------------------
                # Uptake is Arg
                # ----------------------------------------
                uptakes[patient] = patients[patient].uptake(7)
            # ----------------------------------------
            # Create a file with all results
            # ----------------------------------------
            excelfile = Excel()
            groups = list(set([key.split('-')[0] for key in uptakes.keys()]))
            groups.sort()
            patient_data = [['patient'] + patients[uptakes.keys()[0]].samples]
            for group in groups:
                keys = []
                group_data = []
                for key in uptakes.keys():
                    if group in key:
                        keys.append(key)
                        group_data.append(uptakes[key])
                group_data = np.array(group_data)
                mean = group_data.mean(axis = 0)
                std = group_data.std(axis = 0)

                for key in keys:
                    patient_data.append([key] + uptakes[key])
                patient_data.append(['' for _ in xrange(len(mean)+1)])
                patient_data.append(['Mean'] + list(mean))
                patient_data.append(['Std'] + list(std))
                patient_data.append(['' for _ in xrange(len(mean)+1)])
                patient_data.append(['' for _ in xrange(len(mean)+1)])

            patient_data = np.array(patient_data)
            excelfile.add_sheet('patient_data', patient_data.T)
            excelfile.save_excel_file(folder + 'patients.xls')
            # ----------------------------------------
            # Finish
            # ----------------------------------------
            self.button_data.Label = 'Saved!'
        #
        openDirDialog.Destroy()
        #
        return None
    ### ========================================================================
    ### ========================================================================
    ### ========================================================================
