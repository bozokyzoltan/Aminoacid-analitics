# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 15:30:29 2015

@author: zoltan
"""

import os
import numpy as np


class AminoacidAnalyzes_file(object):
    """
    """
    def __init__(self, datafilename):
        """
        """
        self.novalue = -9999
        self.read_datafile(datafilename)
        #
        return None
    ### ==================================================================== ###
    def read_datafile(self, datafilename):
        """
        """
        with open(datafilename, 'r') as datafile:
            self.lines = [line.split('\t') for line in datafile.read().replace('"', '').splitlines()]
        #
        self.filename = os.path.basename(datafilename)
        #
        for index, line in enumerate(self.lines):
            if 'SampleName' in line:
                self.sample_id = line[1]
            if 'Customer_Info' in line:
                self.sample = line[1].split()[-1]
            if 'Amt_Taken' in line:
                self.amount_taken = float(line[1].split()[0])
            if '#' in line:
                start_index = index + 1
        #
        self.patient = '-'.join(self.sample.split('-')[:2])
        self.a_or_b = [0, 1]['B' in self.sample.upper()]
        self.tag = self.sample.split('-')[-1].replace('B','')
        #
        self.aa = []
        self.amount = []
        self.mask = []
        self.conc = []
        #
        for line in self.lines[start_index : start_index + 23]:
            self.aa.append(line[1])
            usefull_data = True
            try:
                self.amount.append(float(line[4]))
            except ValueError:
                self.amount.append(self.novalue)
                usefull_data = False
            except IndexError:
                self.amount.append(self.novalue)
                usefull_data = False
            self.mask.append(usefull_data)
            self.conc.append(self.amount[-1] / self.amount_taken)
        # -------------------
        # Convert to numpy array
        # -------------------
        self.aa = np.array(self.aa)
        self.amount = np.array(self.amount, dtype = np.float32)
        self.conc = np.array(self.conc, dtype = np.float32)
        self.mask = np.array(self.mask, dtype = np.bool8)
        #
        return None
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###

class Patient(object):
    """
    """
    def __init__(self):
        """
        """
        self.data = {}
        self.novalue = -9999
        self.samples = ['t=0', '0', '1', '2', '3', '4', '5']
        #
        return None
    ### ==================================================================== ###
    def add_data(self, aa_object):
        """
        """
        if aa_object.tag not in self.data:
            self.data[aa_object.tag] = []
        self.data[aa_object.tag].append(aa_object)
        #
        return None
    ### ==================================================================== ###
    def calculate_average(self):
        """
        """
        self.aa = self.data[self.data.keys()[0]][0].aa
        self.mean = np.empty((7, len(self.aa)), dtype = np.float32)
        self.std = np.empty((7, len(self.aa)), dtype = np.float32)
        self.mask = np.zeros((7, len(self.aa)), dtype = np.bool8)
        for sample_index in xrange(len(self.samples)):
            if self.samples[sample_index] in self.data:
                mean_temp = []
                mask_temp = []
                for paralell in xrange(len(self.data[self.samples[sample_index]])):
                    mean_temp.append(self.data[self.samples[sample_index]][paralell].conc)
                    mask_temp.append(self.data[self.samples[sample_index]][paralell].mask)
                mean_temp = np.array(mean_temp, dtype = np.float32)
                mask_temp = np.array(mask_temp, dtype = np.bool8)
                for j in xrange(mean_temp.shape[1]):
                    values = mean_temp[:, j][mask_temp[:, j]]
                    if len(values) > 0:
                        self.mask[sample_index, j] = True
                        self.mean[sample_index, j] = mean_temp[:, j][mask_temp[:, j]].mean()
                        self.std[sample_index, j] = mean_temp[:, j][mask_temp[:, j]].std()
                    else:
                        self.mask[sample_index, j] = False
                        self.mean[sample_index, j] = self.novalue
                        self.std[sample_index, j] = self.novalue
        #
        self.normalize()
        self.calculate_uptake_for_each()
        #
        return None
    ### ==================================================================== ###
    def normalize(self):
        """
        """
        self.normed = {}
        for col in xrange(-1, 1):
            self.normed[col] = {}
            #
            mask = self.mask.copy()
            mean = self.mean.copy()
            #
            if col >= 0:
                # Normalize to the first column
                for i in xrange(mask.shape[0]):
                    for j in xrange(mask.shape[1]):
                        if mask[col, j]:
                            mean[i, j] = self.mean[col, j] - mean[i, j]
                        else:
                            mask[i, j] = False
                self.normed[col][-1] = (mask, mean)
            #
            # Normalize to each row
            for row_number in xrange(len(self.aa)):
                new_mask = mask.copy()
                new_mean = mean.copy()
                for i in xrange(new_mask.shape[0]):
                    for j in xrange(new_mask.shape[1]):
                        if ((new_mask[i, row_number]) and (mean[i, row_number] != 0.0)):
                            new_mean[i, j] *= 1.0 / mean[i, row_number]
                        else:
                            new_mask[i, j] = False
                self.normed[col][row_number] = (new_mask, new_mean)
        #
        return new_mean
    ### ==================================================================== ###
    def calculate_uptake_for_each(self):
        """
        """
        self.uptake = []
        t0_mean = 0.0
        t0_mask = False
        for sample_index in xrange(len(self.samples)):
            if self.samples[sample_index] in self.data:
                value_temp = []
                mask_temp = []
                for paralell in xrange(len(self.data[self.samples[sample_index]])):
                    value_temp.append(self.data[self.samples[sample_index]][paralell].conc)
                    mask_temp.append(self.data[self.samples[sample_index]][paralell].mask)
                value_temp = np.array(value_temp, dtype = np.float32)
                mask_temp = np.array(mask_temp, dtype = np.bool8)
                if sample_index == 0:
                    for j in xrange(1, value_temp.shape[1]):
                        values = value_temp[:, j][mask_temp[:, j]]
                        print values
                            if len(values) > 0:
                                t0_mean = values[:, 0].mean()
                                t0_mask = True
    
                    for j in xrange(value_temp.shape[1]):
                        values = value_temp[:, j][mask_temp[:, j]]
                        if len(values) > 0:
                            pass
        #
        exit()
        #
        return None
    ### ==================================================================== ###
    def uptake(self, row):
        """
        """
        data = self.normed[0][-1]
        arg = []
        for i in xrange(len(self.samples)):
            if data[0][i, row]:
                arg.append(data[1][i, row])
            else:
                arg.append('')
        #
        return arg
    ### ==================================================================== ###
    def raw_data(self):
        """
        """
        raw_values = [['Filename', 'Sample_ID', 'Sample'] + list(self.aa) + ['', 'Amount taken /ul']]
        conc_values = [['Filename', 'Sample_ID', 'Sample'] + list(self.aa)]
        amount_taken = ['Amount taken /ul']
        for sample_index in xrange(len(self.samples)):
            if self.samples[sample_index] in self.data:
                for paralell in xrange(len(self.data[self.samples[sample_index]])):
                    mask = self.data[self.samples[sample_index]][paralell].mask
                    filename = self.data[self.samples[sample_index]][paralell].filename
                    sample_id = self.data[self.samples[sample_index]][paralell].sample_id
                    sample = self.samples[sample_index]
                    raw_temp = [filename, sample_id, sample]
                    conc_temp = [filename, sample_id, sample]
                    for i in xrange(len(mask)):
                        if mask[i]:
                            raw_temp.append(self.data[self.samples[sample_index]][paralell].amount[i])
                            conc_temp.append(self.data[self.samples[sample_index]][paralell].conc[i])
                        else:
                            raw_temp.append('')
                            conc_temp.append('')
                    raw_temp.append('')
                    raw_temp.append(self.data[self.samples[sample_index]][paralell].amount_taken)
                    raw_values.append(raw_temp)
                    conc_values.append(conc_temp)
                    amount_taken.append(self.data[self.samples[sample_index]][paralell].amount_taken)
        #
        return conc_values, raw_values, amount_taken
    ### ==================================================================== ###
    def mean_std(self):
        """
        """
        mean = [['Sample', 'Aminoacid'] + list(self.aa)]
        std = [['Sample', 'Aminoacid'] + list(self.aa)]
        for i in xrange(self.mask.shape[0]):
            temp_mean = [self.samples[i], 'Mean']
            temp_std = [self.samples[i], 'Std']
            for j in xrange(self.mask.shape[1]):
                if self.mask[i, j]:
                    temp_mean.append(self.mean[i, j])
                    temp_std.append(self.std[i, j])
                else:
                    temp_mean.append('')
                    temp_std.append('')
            mean.append(temp_mean)
            std.append(temp_std)
        #
        return mean + [''] + std
    ### ==================================================================== ###
    def normalized(self, col, row):
        """
        """
        mask = self.normed[col][row][0]
        mean = self.normed[col][row][1]

        value = [[self.aa[row] + ' normalized', 'Sample'] + list(self.aa) + ['']]
        for i in xrange(mask.shape[0]):
            temp_value = ['Normalized', self.samples[i]]
            for j in xrange(mask.shape[1]):
                if mask[i, j]:
                    temp_value.append(mean[i, j])
                else:
                    temp_value.append('')
            temp_value.append('')
            value.append(temp_value)
        #
        return value
    ### ==================================================================== ###
    ### ==================================================================== ###
    ### ==================================================================== ###
