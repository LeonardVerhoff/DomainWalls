"""
Script to create supercells of LN and LT in orthorhombic cell along x- or y-direction.
Also: Supercell of hexagonal supercell along z-direction
Optional: different domains.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import pymatgen.core as mg
from pymatgen.io.vasp import Poscar
import argparse

"""Global arguments"""
num_neighbors = {'Li': 3, 'Nb': 6, 'Ta': 6} 
nn_dist = 3.0 # radius of sphere where we search for neighbors 
width_scale = 10 # dictating steepness of tanh


def main():

    print('')

    #argument parsing
    parser = argparse.ArgumentParser(description='Multiply given orthorombic unit cell of LN or LT. \
                                     Additionally, 2 domain walls can be placed in the cell \
                                     at 0.25 and 0.75 in either x or y direction.\
                                     DW parallel to x axis is X-type; parallel to y axis is Y-type.', 
                                    epilog='Written by Leonard M. Verhoff 2023\n     Last modified: 03/2024', 
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--infile', '-i', required=True, type=str,
                        help='The POSCAR of unit cell we want to modify.')
    parser.add_argument('--outfile', '-o', default='./supercell.vasp', type=str,
                        help='Where to store (i.e. directoy+name) of produced cell.')
    parser.add_argument('--xmult', '-x', default=1, type=int, 
                        help='Number of multiplaction in x-direction')
    parser.add_argument('--ymult', '-y', default=1, type=int, 
                        help='Number of multiplaction in y-direction')
    parser.add_argument('--zmult', '-z', default=1, type=int, 
                        help='Number of multiplaction in z-direction')
    parser.add_argument('--walltype', '-w', default='None', type=str, choices=['None', 'X', 'Y', 'XY'],
                        help='If we want to create opposite poled domains. X-type: parallel to x axis, \
                            Y-type parallel to y-axis. XY-type: parallel to XY-plane.')
    parser.add_argument('--tanh', '-t', default='False', 
                        help='Decide wether to already assume tanh-shift.')
    parser.add_argument('--xtype', '-X', default='None', choices=['None', 'I', 'II'],
                        help='There are 2 types of X-walls. With classical method, we get only average. Hence we use truncated bulk approach. \n \
                        I: Up -> down -> vacuum; II: Down -> Up -> vacuum \n \
                            Watch out: This does not work with tanh-mode!')
    parser.add_argument('--vacuum', '-v', default=0.0, type=float, 
                        help='Distance of vacuum between 2 blocks in Angstroem.')
    parser.add_argument('--fix', '-f', default=0, type=int, nargs=3, 
                        help='Weather to fix atoms at places (x, y, and z direction). Please provide three numbers (0s and 1s).')
    parser.add_argument('--dfix', '-df', default=0.05, type=float,
                        help='Determines the width in the center of domain that is held constant in units of the c-vector. Only implemented for XY-walls. \
                        (if -f != 0)' )
    parser.add_argument('--fixO', '-fO', default=False, type=int,
                        help='decides weather the whole oxygen sublattice shall be fixed.' )
    

    args = parser.parse_args()

    
    # scale unit vectors by the factors specified in the input
    scale = np.array([args.xmult, args.ymult, args.zmult]) 

    # load structure
    poscar_supercell = Poscar.from_file(args.infile)
    supercell = poscar_supercell.structure

    # create supercell
    supercell.make_supercell(scale)

    # list containing selective dynamics tags for all atoms
    fix_list = len(supercell)*[[1, 1, 1]]
    fix_ctr = 0 # count, how many atoms we fix

    # perfrom only, if we want to create domains
    if args.walltype != 'None':

        # set bound for upper and lower displacements
        lower = 0.25
        upper = 0.75

        # determine axis orthogonal to DW
        if args.walltype == 'Y':
            index = 0

        elif args.walltype == 'XY':
            if args.tanh == True:
                  print("Tanh not implemented for XY. Sorry.")
                  return 0
            #print("Bout to implement for XY.")
            #return 0
            index = 2

        # handeling different X-wall types
        elif args.walltype == 'X':
            index = 1
            if args.xtype== 'II': # type II: invert poling in 1st part
                lower = 0
                upper = 0.5
            elif args.xtype== 'I': # type I: invert poling in 2nd part
                lower = 0.5
                upper = 1.0 
            elif args.xtype != 'None':
                print('please enter valid x-wall type!')
                return 0
        
        else:
            print('Please enter valid walltype!')
            return 0
        
        print(upper, lower)
        
        # iterate over all atoms. If Cation is in central domain, switch its polarization
        for (site_id, site) in enumerate(supercell.sites):
            species = site.specie.symbol

            if species in num_neighbors.keys():
            
                ########## Two domains with poling up and down (step function) ##########
                if args.tanh == 'False':
                    if site.frac_coords[index] >= lower and site.frac_coords[index] < upper:
                        num_nn = num_neighbors[species]
                        # get all neighbors in a sphere of radius nn_dist
                        all_neighbors = supercell.get_neighbors( site, nn_dist )
                        if len(all_neighbors) < 6:
                                print(f'error with {species}! Only found {len(all_neighbors)} neighbors')
                        # determine nearest neighbors. 3 for Li, 6 for Nb/Ta
                        mask = np.argsort([s2.distance(site) for s2 in all_neighbors])
                        # reference: oxygen layer for Li, cage center for Nb/Ta
                        reference = np.mean( [all_neighbors[i].frac_coords[2] for i in mask[0:num_nn]] )
                        displ = site.frac_coords[2] - reference
                        site.frac_coords[2] -= 2.0*displ 

                    

                ########### Two domains with poling up and down assuming tanh! (-> faster convergence?) ##########
                elif args.tanh == 'True':
                        #print(args.tanh)
                    #if site.frac_coords[index] >= 0.25 and site.frac_coords[index] < 0.75:
                        num_nn = num_neighbors[species]
                        # get all neighbors in a sphere of radius nn_dist
                        all_neighbors = supercell.get_neighbors( site, nn_dist )
                        if len(all_neighbors) < 6:
                                print(f'error with {species}! Only found {len(all_neighbors)} neighbors')
                        # determine nearest neighbors. 3 for Li, 6 for Nb/Ta
                        mask = np.argsort([s2.distance(site) for s2 in all_neighbors])
                        # reference: oxygen layer for Li, cage center for Nb/Ta
                        reference = np.mean( [all_neighbors[i].frac_coords[2] for i in mask[0:num_nn]] )
                        p0 = site.frac_coords[2] - reference
                        fac = width_scale#/(args.xmult*args.ymult)
                        if site.frac_coords[index] <= 0.5:
                            site.frac_coords[2] = reference - p0 * np.tanh( fac*(site.frac_coords[index] - 0.25) )

                        elif site.frac_coords[index] > 0.5:
                            site.frac_coords[2] = reference + p0 * np.tanh( fac*(site.frac_coords[index] - 0.75) )
                        
                        else:
                             print('Error')
                             return 0
                        
            # if species is not oxygen, we may include it in fix_list:
            if args.fix != 0:
                if args.walltype == 'XY':
                        # in case the reference z-pos is in (center-df/2, center+df/2) -> add index to fix_list!
                    z_pos = site.frac_coords[2]
                    if (z_pos < args.dfix/2 and z_pos > -args.dfix/2) or (z_pos < 0.5+args.dfix/2 and z_pos > 0.5-args.dfix/2): 
                        # for first domain or 2nd domain
                        fix_list[site_id] = [1-args.fix[0], 1-args.fix[1], 1-args.fix[2]]
                        fix_ctr += 1    
                    

                else:
                        print('Sorry, fixing is only implemented for XY-walls (orthogonal to z direction)')

            if species == 'O' and args.fixO == 1:
                fix_list[site_id] = [0, 0, 0] 
                        
         

    supercell_ = Poscar(supercell)  

    # set selective dynamics tags from list
    if args.fix != 0:
        if args.walltype == 'XY': 
            supercell_.selective_dynamics = fix_list
            print(f'Fixed {fix_ctr} cations!')
            if args.fixO == 1:
                print('fixed all oxygen anions.')
             
             
    supercell = supercell_.structure
    
    # if we have special X-treatment, we need to create vacuum region:
    if args.walltype == 'X' and (args.xtype == 'I' or args.xtype == 'II'):
         # Get the current lattice parameters
        lattice = supercell.lattice
        current_lengths = lattice.lengths
        current_angles = lattice.angles

        # Create the modified lattice with added vacuum region along y-direction
        modified_lattice = mg.Lattice.from_parameters(current_lengths[0], current_lengths[1]+args.vacuum, current_lengths[2], *current_angles)
        # Create new supercell with modified lattice
        supercell = mg.Structure(modified_lattice, supercell.species, supercell.cart_coords, coords_are_cartesian=True)
         

    poscar = Poscar(supercell)
    poscar.write_file(args.outfile)               

    print(f'------------------------------\nWrote new {args.xmult}x{args.ymult}x{args.zmult} cell to {args.outfile}.\n------------------------------\n')


if __name__ == "__main__":
    main()