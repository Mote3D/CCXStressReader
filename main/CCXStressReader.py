#!/usr/bin/env python
# -*- coding: utf-8 -*-

__copyright__ = "Copyright (C) 2021 Henning Richter"
__email__ = "mote3d@quantentunnel.de"
__version__ = "1.0"
__license__ = """
    This library is free software; you can redistribute it and/or 
    modify it under the terms of the GNU Lesser General Public 
    License as published by the Free Software Foundation; either 
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    along with this library; if not, see <http://www.gnu.org/licenses/>.
"""

"""
CCXStressReader reads element variable output such as stresses and strains 
from .dat files created by the FOSS finite-element software CalculiX [G. Dhondt, 
K. Wittig, http://www.calculix.de/] or the FreeCAD [J. Riegel, W. Mayer, 
Y. van Havre, http://www.freecadweb.org/] FEM module and computes the minimum, 
maximum and arithmetic mean of Mises equivalent stress, total effective strain 
and equivalent plastic strain. The results are written to a .txt file.

The .dat file stores element variable output at the elements' integration 
points, which generally is more accurate than output extrapolated to the nodes. 
Evaluating element variable output at the integration points is therefore 
preferred, particularly in stress analyses with non-linear material behaviour.

Element variable output to the .dat file can be activated by adding the 
following command to the CalculiX .inp file:

*EL PRINT, ELSET=Eall
S, E, PEEQ

CCXStressReader has been tested with Python 3.7.2 and NumPy 1.16.2.
"""


import numpy as np
import sys, time


def read_input(fname):
    """Read element variable output at integration points from .dat file"""
    with open(fname, 'r') as f:
        lines = [line.strip().split() for line in f if line.strip()]

    keylist = [' '.join(cline[0:3]) for cline in lines if not cline[0].isdigit()]
    valuelist = [' '.join(cline[0:3]) for cline in lines]
    otherlist = [cline for cline in keylist if cline not in ['stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz)',
                 'strains (elem, integ.pnt.,exx,eyy,ezz,exy,exz,eyz)', 'equivalent plastic strain']]

    # Find indices:
    if 'stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz)' in keylist:
        Sindex = valuelist.index('stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz)')
    else:
        print("\nStress output S not found.\n")
        time.sleep(1)
        sys.exit()

    if 'strains (elem, integ.pnt.,exx,eyy,ezz,exy,exz,eyz)' in keylist:
        Eindex = valuelist.index('strains (elem, integ.pnt.,exx,eyy,ezz,exy,exz,eyz)')
    else:
        Eindex = np.nan
        print("\nStrain output E not found.\n")

    if 'equivalent plastic strain' in keylist:
        PEEQindex = valuelist.index('equivalent plastic strain')
    else:
        PEEQindex = np.nan
        print("\nEquivalent plastic strain output PEEQ not found.\n")

    if otherlist:
        Oindexlist = []
        for i in range(0, len(otherlist), 1):
            Oindexlist.append(valuelist.index(otherlist[i]))
        Oindex = np.nanmin(Oindexlist)
    else:
        Oindex = np.nan

    indexarr = np.array([Sindex, Eindex, PEEQindex, Oindex, np.nan])
    indexarrsort = np.argsort(indexarr)

    # Read element variable output:
    if np.isnan(indexarr[indexarrsort[np.where(indexarrsort==0)[0][0]+1]]):
        Sarr = np.array(lines[(Sindex+1):len(lines)], dtype=float)
    else:
        Sarr = np.array(lines[(Sindex+1):int(indexarr[indexarrsort[np.where(indexarrsort==0)[0][0]+1]])], dtype=float)

    if np.isnan(Eindex):
        Earr = np.full((Sarr.shape[0], 3), np.nan)
    elif np.isnan(indexarr[indexarrsort[np.where(indexarrsort==1)[0][0]+1]]):
        Earr = np.array(lines[(Eindex+1):len(lines)], dtype=float)
    else:
        Earr = np.array(lines[(Eindex+1):int(indexarr[indexarrsort[np.where(indexarrsort==1)[0][0]+1]])], dtype=float)

    if np.isnan(PEEQindex):
        PEEQarr = np.full((Sarr.shape[0], 3), np.nan)
    elif np.isnan(indexarr[indexarrsort[np.where(indexarrsort==2)[0][0]+1]]):
        PEEQarr = np.array(lines[(PEEQindex+1):len(lines)], dtype=float)
    else:
        PEEQarr = np.array(lines[(PEEQindex+1):int(indexarr[indexarrsort[np.where(indexarrsort==2)[0][0]+1]])], dtype=float)
    return (Sarr, Earr, PEEQarr)


def compute_eqstress(S):
    """Compute Mises equivalent stress"""
    MISESarr = np.column_stack((S[:,0], S[:,1], np.sqrt(0.5*((S[:,2]-S[:,3])*(S[:,2]-S[:,3])
                              +(S[:,3]-S[:,4])*(S[:,3]-S[:,4])+(S[:,4]-S[:,2])*(S[:,4]-S[:,2]))
                              +3.0*(S[:,5]*S[:,5]+S[:,6]*S[:,6]+S[:,7]*S[:,7]))))
    return MISESarr


def compute_eqstrain(E):
    """Compute total effective strain"""
    EEQarr = np.column_stack((E[:,0], E[:,1], (2.0/3.0)*np.sqrt(0.5*((E[:,2]-E[:,3])*(E[:,2]-E[:,3])
                            +(E[:,3]-E[:,4])*(E[:,3]-E[:,4])+(E[:,4]-E[:,2])*(E[:,4]-E[:,2]))
                            +3.0*(E[:,5]*E[:,5]+E[:,6]*E[:,6]+E[:,7]*E[:,7]))))
    return EEQarr


def write_txtfile(fname, RESarr):
    """Write results to .txt file"""
    MINarr = np.amin(RESarr, axis=0)
    MAXarr = np.amax(RESarr, axis=0)
    MEANarr = np.mean(RESarr, axis=0)

    with open(fname[:-4]+'_IntPtOutput.txt', 'w') as d:
        fmt = '{}{}{}'.format('%12i', '%9i', ' '.join(['%16.4e']*3))
        hdr = '{0}Elem.    Int.Pt.{1}MISES{2}EEQ{3}PEEQ'.format(' '*5, ' '*9, ' '*14, ' '*13)
        np.savetxt(d, RESarr, fmt=fmt, delimiter=' ', newline='\n', header=hdr)
    
        d.write('\n')
        d.write('{0}Minimum{1:25.4e}{2:17.4e}{3:17.4e}\n'.format(' '*5, MINarr[2], MINarr[3], MINarr[4]))
        d.write('{0}Maximum{1:25.4e}{2:17.4e}{3:17.4e}\n'.format(' '*5, MAXarr[2], MAXarr[3], MAXarr[4]))
        d.write('{0}Mean (arith.){1:19.4e}{2:17.4e}{3:17.4e}\n'.format(' '*5, MEANarr[2], MEANarr[3], MEANarr[4]))
        d.write('\n')
    print("\nResults successfully written to file '{}'.\n".format(fname[:-4]+'_IntPtOutput.txt'))


def main():

    # Specify path to .dat file:
    fname = input("\nPlease enter path and name of .dat file:\n")

    if not fname[-4:] == ".dat":
        print("\nFile is invalid, please select .dat file.\n")
        time.sleep(1)
        sys.exit()

    # Read element variable output from .dat file:
    (S, E, PEEQ) = read_input(fname)

    # Compute Mises equivalent stress:
    MISES = compute_eqstress(S)

    # Compute total effective strain:
    if not np.isnan(np.sum(E[:,2], axis=0)):
        EEQ = compute_eqstrain(E)
    else:
        pass

    # Format output:
    RESarr = np.column_stack((MISES, EEQ[:,2], PEEQ[:,2]))

    # Write results to .txt file:
    write_txtfile(fname, RESarr)


if __name__ == "__main__":
    main()
