# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 13:09:11 2012

@author: balarsen
"""

def listToEnumerate(inList, startingNum=1):
    """
    starting from a python list return a string that is the RST equlivant of
    the list in enumerated list
    """
    outVal = ''
    for i, val in enumerate(inList):
        outVal += '{0}. {1}\n'.format(i+startingNum, val)
    outVal += '\n'
    return outVal

def listToList(inList):
    """
    starting from a python list return a string that is the RST equlivant of
    the list in a bulleted list
    """
    outVal = ''
    for i, val in enumerate(inList):
        outVal += '- {0}\n'.format(val)
    outVal += '\n'
    return outVal

def listToTable(data, header='', title=''):
    """
    starting from a python list return a string that is the RST equlivant of
    the list in rst table format

    based loosly on
    http://stackoverflow.com/questions/11347505/what-are-some-approaches-to-outputting-a-python-data-structure-to-restructuredte
    see
    http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html#tables
    """
    numcolumns = len(data[0])
    colsizes = [str(max(len(str(r[i])) for r in data)) for i in range(numcolumns)]
    code_header = '.. csv-table:: {0}'.format(title)
    header = '\t:header: ' + ', '.join(header)
    widths = '\t:widths: ' + ', '.join(colsizes)

    data_out = ''
    for row in data:
        data_out += '\t' + ', '.join([str(v) for v in row]) + '\n'
    output = code_header + '\n' + header + '\n' + widths + '\n\n' + data_out
    return output







