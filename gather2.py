#    gather2.py : Data Output Support for COVID-19 Simulation
#    Copyright (C) 2020 Philip T. Gressman <gresssman@math.upenn.edu> and Jennifer R. Peck <jpeck1@swarthmore.edu>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.



import os
import re
import math
import sys,time

def custom_sum(item1,item2):
    if item1 is None:
        return item2
    if type(item1) == list:
        if len(item1) > 0 and (type(item2) != list or len(item1) != len(item2)):
            print(item1,item2)
            print(len(item1),len(item2))
            raise Exception('not lists of matching size')
        length = max(len(item1),len(item2))
        result = [0] * length
        for index in range(length):
            if index < len(item1):
                result[index] += item1[index]
            if index < len(item2):
                result[index] += item2[index]
        return result
    return item1 + item2

def decompress(shortlist,length):
    result = []
    if len(shortlist) == 0:
        return [0] * length
    for item in shortlist:
        if type(item) == list:
            result += [item[0]] * item[1]
        else:
            result.append(item)
    return result
def compare(dict1,dict2):
    for key in dict1:
        if key not in dict2:
            return False
        elif dict1[key] != dict2[key]:
            return False
    for key in dict2:
        if key not in dict1:
            return False
    return True

class Summarizer(object):
    def __init__(self,lists_of_numbers):
        self.column_titles = []
        self.row_titles = []
        self.lists_of_numbers = lists_of_numbers
        self.summary = []
    def quantiles(self,quantilelist):
        self.summary = summarize(self.lists_of_numbers)

def summarize(listthings,quantiles):
    if type(listthings[0]) == list:
        result = []
        for index in range(len(listthings[0])):
            itemsatspot = []
            for preind in range(len(listthings)):
                itemsatspot.append(listthings[preind][index])
            result.append(summarize(itemsatspot,quantiles))
        return result
    mylist = listthings
    result = []
    mylist.sort()
    for quant in quantiles:
        frac_index = (len(mylist) - 1) * quant
        lower = math.floor(frac_index)
        upper = math.ceil(frac_index)
        if lower == upper:
            result.append(mylist[lower])
        else:
            vlow = mylist[lower]
            vhigh = mylist[upper]
            interp = vlow * (upper - frac_index) + vhigh * (frac_index - lower)
            result.append(interp)
    return result

class DataCollector(object):
    def __init__(self,where):
        self.data_found = {}
        self.data_groups = []
        self.keys = {}
        self.reports = []
        self.quantiles = [0.95,0.75,0.5,0.25,0.05]
        self.generated_report = None
        self.keep_raw = False
        self.gather_reports(where)
        #self.pool_data()
        #self.load_all()
        self.populate_keys()
    def register_report(self,name,pattern,function=None):
        self.reports.append([name,pattern,function])
    def gather_reports(self,runobject):
        self.data_found = {}
        self.parameters_found = []
        self.data_groups = []
        found = False
        structure = runobject.recorder.compress()
        fullname = 'current_run'
        if type(structure) == dict and '_information' in structure:
            length = structure['_information']['run_days'] + 1
            for key in structure:
                if key != '_information':
                    for index,listdata in enumerate(structure[key]):
                        structure[key][index] = decompress(listdata,length)
            if 'filename' in structure['_information']:
                del structure['_information']['filename']
            self.data_found[fullname] = structure
            print('===== Data Ready for Analysis.')
            found = True

        if found:
            pooled = False
            for index in range(len(self.parameters_found)):
                if compare(self.parameters_found[index],structure['_information']):
                    self.data_groups[index].append(fullname)
                    pooled = True
            if not pooled:
                self.parameters_found.append(structure['_information'])
                self.data_groups.append([fullname])

    def populate_keys(self):
        self.keys = {}
        for item in self.data_found:
            for subitem in self.data_found[item]:
                if subitem[0] != '_':
                    self.keys[subitem] = True
    def matching_keys(self,mydict):
        if type(mydict) == list:
            for item in mydict:
                if item not in self.keys:
                    print('!!!!! WARNING: ' + item + ' not known')
            return mydict
        found_matches = []
        match_profile = {}
        for key in mydict:
            match_profile[key] = [False,False]
        for key in self.keys:
            matches = True
            if type(key) == tuple:
                for item,value in mydict.items():
                    if (item in key and not value) or (item not in key and value):
                        matches = False
                        match_profile[item][0] = True
                    else:
                        match_profile[item][1] = True
            else:
                matches = False
            if matches:
                found_matches.append(key)
        for key in match_profile:
            if not match_profile[key][0] and not match_profile[key][1]:
                print('!!!!! WARNING: No compound keys found here.')
                print(key)
                print(self.keys)
            elif not match_profile[key][0]:
                print('!!!!! WARNING: All keys match ' + key + ' : ' + str(mydict[key]))
            elif not match_profile[key][1]:
                print('!!!!! WARNING: No keys match ' + key + ' : ' + str(mydict[key]))

        return found_matches
    def collect(self,name,matching_keys,apply_function=None):
        final_answer = []
        for grouping in self.data_groups:  # Run is a list of identically configured runs
            medium_answer = []
            count = 0
            scname = ''
            result = None
            sumlists = []
            offset = 0
            for runno in grouping:
                if 'scenario_name' in self.data_found[runno]['_information']:
                    scname = self.data_found[runno]['_information']['scenario_name']
                for key in matching_keys:
                    if key in self.data_found[runno]:
                        for index,timeseries in enumerate(self.data_found[runno][key]):
                            if len(sumlists) <= offset+index:
                                sumlists.append([])
                            sumlists[offset+index] = custom_sum(sumlists[offset+index],timeseries)
                offset = len(sumlists)
            count += len(sumlists)
            if len(sumlists) == 0:
                sumlists = [[0]]
            for index in range(len(sumlists)):
                if apply_function is not None:
                    sumlists[index] = apply_function(sumlists[index])
            medium_answer = sumlists
            if len(medium_answer) > 0:
                remember = medium_answer
                medium_answer = summarize(medium_answer,self.quantiles)
                result_dict = {'scenario_name' : scname, 'run_pools' : grouping, 'count' : count, 'quantiles' : medium_answer}
                if self.keep_raw:
                    result_dict['raw'] = remember
                final_answer.append(result_dict)
                print('===== Quantity:',name,'Scenario:',scname,'Count:',count)
        return final_answer
    def generate_reports(self):
        final_report = {}
        for name,pattern,function in self.reports:
            keyset = self.matching_keys(pattern)
            result = self.collect(name,keyset,function)
            final_report[name] = result
        return final_report
    def generate_csv(self):
        csv = CSVFile()
        result = self.generate_reports()
        self.generated_report = result
        column = 0
        for key in result:
            csv.set(0,column,'scenario')
            column += 1
            for index,quantile in enumerate(self.quantiles):
                csv.set(0,column,key + ' ' + str(quantile))
                row = 1
                for item in result[key]:
                    data = item['quantiles']
                    if len(data) > 0:
                        if type(data[0]) == list:
                            for line in data:
                                csv.set(row,column,line[index])
                                if index == 0:
                                    csv.set(row,column-1,item['scenario_name'])
                                row += 1

                        else:
                            csv.set(row,column,data[index])
                            if index == 0:
                                csv.set(row,column-1,item['scenario_name'])
                            row += 1
                column += 1
        return csv


class CSVFile(object):
    def __init__(self):
        self.cells = {}
        self.rowhome = 0
        self.colhome = 0
        self.col_lock = -1
        self.atrow = 0
        self.atcol = 0
        self.maxrow = 0
        self.maxcol = 0
    def rehome(self):
        self.rowhome = 0
        self.colhome = 0
        self.col_lock = -1
    def column_length(self,colnumber):
        howhigh = 0
        for row in self.cells:
            if colnumber in self.cells[row]:
                howhigh = row
        return howhigh + 1
    def read(self,row,column):
        return self.cells[row][column]
    def output(self):
        result = ''
        for row in range(self.maxrow+1):
            if row >= self.rowhome:
                thisline = ''
                for col in range(self.maxcol+1):
                    if col >= self.colhome and (self.col_lock == -1 or col < self.colhome + self.col_lock):
                        contents = self.get(row,col)
                        if contents is None:
                            contents = ''
                        elif type(contents) == str:
                            contents = '"' + contents + '"'
                        else:
                            contents = str(contents)
                        thisline += contents + ','
                thisline = thisline[0:len(thisline)-1] + '\n'
            result += thisline
        if len(result) > 0:
            result = result[0:len(result)-1]
        return result

    def set_lock(self,row,col,colnumb):
        self.rowhome = row
        self.colhome = col
        self.collock = colnumb
        self.atrow = row
        self.atcol = col
    def next(self):
        self.atcol += 1
        if self.atcol == self.colhome + self.col_lock:
            self.atcol = self.colhome
            self.atrow += 1
    def next_row(self):
        self.atrow = self.atrow + 1
        self.atcol = self.colhome
    def next_column(self):
        self.atcol += 1
    def set(self,row,col,val):
        if row is None:
            row = self.atrow
        if col is None:
            col = self.atcol
        if row not in self.cells:
            self.cells[row] = {}
        self.maxrow = max(row,self.maxrow)
        self.maxcol = max(col,self.maxcol)
        self.cells[row][col]  = val
    def get(self,row=None,col=None):
        if row is None:
            row = self.atrow
        if col is None:
            col = self.atcol
        if row not in self.cells:
            return None
        if col not in self.cells[row]:
            return None
        return self.cells[row][col]
