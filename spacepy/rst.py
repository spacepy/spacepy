# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 13:09:11 2012

@author: balarsen
"""

def listToEnumerate(inList, startingNum=1, makeBlock=True):
    """
    starting from a python list return a string that is the RST equlivant of
    the list in enumerated list
    makeBlock : make the text into a text block in rst if there are \n in the string
    """
    outVal = ''
    for i, val in enumerate(inList):
        if makeBlock:
            val = val.replace('\n', '\n    ') # makes it a text block
        outVal += '{0}. {1}\n'.format(i+startingNum, val)
    outVal += '\n'
    return outVal

def listToList(inList, makeBlock=True):
    """
    starting from a python list return a string that is the RST equlivant of
    the list in a bulleted list
    """
    outVal = ''
    for i, val in enumerate(inList):
        if makeBlock:
            val = val.replace('\n', '\n    ') # makes it a text block
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

def strToHeading(string, level=0):
    """
    return a rst heading from the given string (soplit on spaces)
    document title: === aboce and below (-2)
    document subtitle: --- above and below (-1)
    section 1: ==== (0)
    subsection 1.1: ---- (1)
    subsubsection 1.1.1: ~~~~ (2)
    """
    lens = [len(v) for v in string.split()]
    if level in [0, -2]:
        mark = '='
    elif level in [1, -1]:
        mark = '-'
    elif level == 2:
        mark = '~'
    else:
        raise(ValueError('Bad level given'))
    marks = mark.join([mark*v for v in lens])
    if level < 0:
        out_string = marks + '\n' + string
    else:
        out_string = string
    out_string += '\n' + marks + '\n'
    return out_string






