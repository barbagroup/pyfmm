"""
Filename: simlogger.py
Author: Felipe Cruz V.
Actual Version: 0.001
"""


from numpy import *
import os
import timeimport getopt
import csv


## General constants
POINTS_DIST_SCATTER = 'scatter'
POINTS_DIST_LATTICE = 'lattice'




class CsvLog:
    '''
    Stores information in comma separated value format.
    '''
    
    
    
    def __init__(self, log_name, save_dir):
        ''' Create the file and prepare the file to receive data.'''
        # full file name
        full_path = save_dir + '/' + log_name + '.csv'
        ## Open file in append mode
        self.writer = csv.writer(open(full_path, 'a'))
        self.row = []
        self.rowSize = 0
    
    def addElement(self, element):
        self.row.append(element)
        self.rowSize += 1
    
    
    
    def flushRow(self):
        ''' Save the current working row and sets a new working row '''
        if self.rowSize > 0:
            self.writer.writerow(self.row)
            self.row = []
            self.rowSize = 0
    
    
    
    def addRow(self, row):
        ''' Add a row to the file.'''
        self.writer.writerow(row)
    
    
        

class SimLogger:
    '''
    Logs the simulation results
    '''
    
    def __init__(self, simname, sim_save_dir):
        self.name = simname
        self.directory = sim_save_dir
        
        # Create sim dir if it doesn't exist
        if not os.path.isdir(sim_save_dir):
            os.mkdir(sim_save_dir)
            
    
    
    def csvOutputLog(self, log_name):
        ''' 
        Create a CSV file ready to save data inside the
        simulation logger path.
        '''
        return CsvLog(log_name, self.directory)
    
    
        
    def saveData(self, file_name, dataType, data):
        ''' Save raw data in a file. '''
        if dataType == POINTS_DIST_SCATTER:
            ''' save data as three columns: X, Y, Weight'''
            dataFile = CsvLog(file_name, self.directory)
            Xdata = data[0]                             # X position
            Ydata = data[1]                             # Y position
            Wdata = data[2]                             # weight data
            # save data into a file of 3 columns
            for i in range(len(data[0])):
                row = [Xdata[i], Ydata[i], Wdata[i]]
                dataFile.addRow(row)
    
    
    
