# Script to convert ffnet .net files to .mat files for calculating LANLstar
# Based on gist contributed by Aaron Hendry

# This script assumes that you run in the directory with all the old .net files
# and have gunzipped them.

import glob
import pickle
import numpy as np
import spacepy.datamodel as dm

infiles = glob.glob('*.net')
for infile in infiles:
    # Load the .net file as a pickle
    with open(infile, 'rb') as fh:
        net = pickle.load(fh, encoding='bytes')
    
    # Get the number of input variables and hidden variables
    num_in = len(net.__dict__[b'inno']) + 1
    num_hid = len(net.__dict__[b'hidno'])
    
    # We need to calculate the number and size of the hidden layers. I believe that each model has 
    # either one or two hidden layers. If we have only one layer, it has 20 variables. If we have two
    # layers, the second has 20, and the first makes up the difference
    num_first = num_hid - 20 # Number of variables in the first layer (or zero, if only one layer)
    N = 0
    M = 0

    if num_first > 0:
        two_layers = True
        # The first N weights define the map between the input and the first layer
        N = num_in*num_first
        # In case the weights aren't in numerical order we re-order them:
        ihmap = np.argsort(net.__dict__[b'hidno'][:num_first])
        # Generate the input->hidden weights
        ihweights = net.__dict__[b'weights'][:N].reshape(num_first, num_in)[ihmap, :]
        # The second layer depends on all of the variables in the first layer AND the bias, so we
        # have M weights:   
        M = (num_first+1)*20
        # These also aren't always in order
        hhmap = np.argsort(net.__dict__[b'conec'][N:(N+M), 0].reshape(20, num_first+1)[0, :])
        hhweights = net.__dict__[b'weights'][N:(N+M)].reshape(20, num_first+1)[:, hhmap]
    else:
        two_layers = False
        # The first N weights define the map between the input and the only layer
        N = num_in*20
        # Again, the weights aren't always in numerical order: so re-order them:
        ihmap = np.argsort(net.__dict__[b'hidno'])
        # Generate the input->hidden weights
        ihweights = net.__dict__[b'weights'][:N].reshape(20, num_in)[ihmap, :]
        hhweights = np.zeros((0,21))
    
    # The output depends on the final layer (with 20 variables) and the first input:
    homap = np.argsort(net.__dict__[b'conec'][(N+M):, 0])
    howeights = net.__dict__[b'weights'][(N+M):]
    howeights = howeights[homap]
    
    mfile = {'two_layers': dm.dmarray(np.atleast_1d(two_layers)),
             'inweights': dm.dmarray(net.__dict__[b'eni'][:, 0]),
             'inbias': dm.dmarray(net.__dict__[b'eni'][:, 1]),
             'ihweights': dm.dmarray(ihweights[:, 1:]),
             'ihbias': dm.dmarray(ihweights[:, 0]),
             'hhweights': dm.dmarray(hhweights[:, 1:]),
             'hhbias': dm.dmarray(hhweights[:, 0]),
             'howeights': dm.dmarray(howeights[1:]),
             'hobias': dm.dmarray(np.atleast_1d(howeights[0])),
             'outweight': dm.dmarray(np.atleast_1d(net.__dict__[b'deo'][0][0])),
             'outbias': dm.dmarray(np.atleast_1d(net.__dict__[b'deo'][0][1])),
             }
    mfile = dm.SpaceData(mfile)
    mfile.attrs['Notes'] = 'Converted from {}'.format(infile)
    
    outfile = '.'.join([infile.split('.')[0], 'txt'])
    print('Converting {} -> {}'.format(infile, outfile))
    mfile.toJSONheadedASCII(outfile)