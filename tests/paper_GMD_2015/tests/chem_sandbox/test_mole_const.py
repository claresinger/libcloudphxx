import numpy   as np
import h5py    as h5
import Gnuplot as gp
import math    as mt
import pprint  as pp

import sys
sys.path.insert(0, "../../../../build/bindings/python/")
from libcloudphxx import common as cm

# TODO
# this should all be done as a sum of bins and not an average value
# there is no guarantee that for average value it will be ok
# for parcel model it is checked for sd_conc=1 and it works
# but for here to have it more accurate it should be summed over all bins
# and multiplying b the water mass in each bin and dividing by total mass

# open hdf5 files with data
h5f_ini = h5.File('../../build/tests/chem_sandbox/out_hall_pinsky_stratocumulus/timestep0000000000.h5', 'r')
h5f_spn = h5.File('../../build/tests/chem_sandbox/out_hall_pinsky_stratocumulus/timestep0000010000.h5', 'r')
h5f_end = h5.File('../../build/tests/chem_sandbox/out_hall_pinsky_stratocumulus/timestep0000011800.h5', 'r')

# helper dict for chem names and molar mass
             #name   gas molar mass   aqueous molar mass    label in hdf5     ini    spn   end
help_dict = {
            'H2O2' : [cm.M_H2O2,      cm.M_H2O2,            'H2O2',            0,     0,    0],
            'O3'   : [cm.M_O3,        cm.M_O3,              'O3',              0,     0,    0],
            'SO2'  : [cm.M_SO2,       cm.M_SO2_H2O,         'S_IV',            0,     0,    0],
            'CO2'  : [cm.M_CO2,       cm.M_CO2_H2O,         'C_IV',            0,     0,    0],
            'NH3'  : [cm.M_NH3,       cm.M_NH3_H2O,         'N_III',           0,     0,    0],
            'HNO3' : [cm.M_HNO3,      cm.M_HNO3,            'N_V',             0,     0,    0],
            'H2SO4': [0,              cm.M_H2SO4,           'S_VI',            0,     0,    0]
            }

for key, val in help_dict.iteritems():
 
    name1 = key + "g"

    if key == 'H2SO4':
        name2 = 'chem_S_VI_aq'

        # moles/ug of dry air
        ini = (h5f_ini[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
        spn = (h5f_spn[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
        end = (h5f_end[name2][:] / val[1]).sum() / 76. / 76. * 1e9 

    else:
        name2 = "chem_" + val[2] + "_aq"

        # moles/ug of dry air
        ini = (h5f_ini[name1][:] / val[0] + h5f_ini[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
        spn = (h5f_spn[name1][:] / val[0] + h5f_spn[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
        end = (h5f_end[name1][:] / val[0] + h5f_end[name2][:] / val[1]).sum() / 76. / 76. * 1e9 

    val[3] = ini
    val[4] = spn
    val[5] = end


print "-------------- init vs end of spin-up ----------------------"

for key in ['CO2', 'NH3', 'HNO3', 'SO2', 'O3', 'H2O2', 'H2SO4']:
    
    relative_error = abs(help_dict[key][4] - help_dict[key][3]) / help_dict[key][3]
    print key , " relative error ", relative_error * 100, " %"

print "-------------- init vs end  --------------------------------"

for key in ['CO2', 'NH3', 'HNO3']:
    
    relative_error = abs(help_dict[key][5] - help_dict[key][3]) / help_dict[key][3]
    print key , " relative error ", relative_error * 100, " %"

print "-------------- test react  --------------------------------"

depleted_O3   = help_dict['O3'][3]    - help_dict['O3'][5] 
depleted_H2O2 = help_dict['H2O2'][3]  - help_dict['H2O2'][5] 
gained_S6     = help_dict['H2SO4'][5] - help_dict['H2SO4'][3] 

relative_error = (depleted_O3 + depleted_H2O2 - gained_S6) / gained_S6 

print "created H2SO4 relative error ", relative_error * 100, " %"
