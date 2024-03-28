"""
Script to extract polarization from relaxed structure and write it in a file.

Idea:
    Sample Supercell with a grid, depending on the supercell size. Each gridpoint has a defined, localized dipole moment (displacement*nominal charge).
    Average dipole moment over 2 spatial direction (parallel to DW) and along one sampling point orthogonal to the wall.
    This gives polariaztion depending on distance to wall on the pre-defined grid.
"""

import sys
import numpy as np
import pymatgen.core as mg
from pymatgen.io.vasp import Poscar
import argparse

"""Global arguments"""
atoms_dict = {'Li': [1, 3, 0], 'Nb': [5, 6, 1], 'Ta': [5, 6, 1]} # [charge, number of nn, index] -> index: column for output in pol_grid
nn_dist = 2.7 # use this as radius of sphere to search for neighbors


def main():

    print('')

    #argument parsing
    parser = argparse.ArgumentParser(description='Extract the average polarization along one axis (x, y, or z) \
                        as a vector from CONTCAR file for LN or LT. Output file:\n \
                        distance p_li,x p_li,y, p_li,z p_nb,x p_nb,y, p_nb,z', 
                        epilog='Written by Leonard M. Verhoff 2023\n     Last modified: 03/2024', 
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--infile', '-i', default='./CONTCAR', type=str,
                        help='The CONTCAR file of relaxed mixed-domain structure.')
    parser.add_argument('--outfile', '-o', default='./pol', type=str,
                        help='Where to store (i.e. directoy+name) the extracted list of polarization.')
    parser.add_argument('--xmult', '-x', default=1, type=int, 
                        help='Number of multipliction in x-direction of CONTCAR')
    parser.add_argument('--ymult', '-y', default=1, type=int, 
                        help='Number of multipliction in y-direction of CONTCAR')
    parser.add_argument('--zmult', '-z', default=1, type=int, 
                        help='Number of multipliction in z-direction of CONTCAR')
    parser.add_argument('--system', '-s', default='Nb', type=str,
                        help='Choose wether we have LN (Nb) or LT (Ta)')
    parser.add_argument('--cell', '-c', default='o', type=str,
                        help='Choose the cell symmetry: orthorhombic (o) or hexagonal (h)')
    parser.add_argument('--normalaxis', '-n', default=0, type=int,
                        help='The axis normal to the wall. 0: x-axis (Y-wall), 1: y-axis (X-wall), 2: z-axis (CDW)')
    parser.add_argument('--plot_pz', '-p', default=0, type=int,
                        help='Decide weather to plot Pz for testing purposes.')
    args = parser.parse_args()


    if args.system not in atoms_dict.keys():
          print('Please enter valid system via -s option. Nb for LN, Ta for LT.')
          return 0


    poscar = Poscar.from_file(args.infile)
    structure = poscar.structure

    dimensions = np.array([args.xmult, args.ymult, args.zmult]) # dimensions of supercell in x, y, and z direction 

    # distinguish between orthorhombic (180 DW) or hexagonal cell (CDW):
    if args.cell=='o':
        num_direction = np.array([2, 6, 6]) # in orthorhombic unit cell: 2 lattice points in x direction, 6 lattice points in y direction, 6 in z-direction 
    elif args.cell=='h':
        num_direction = np.array([3, 3, 6]) # in hexagonal unit cell: 3 lattice points in 1st direction, 3 lattice points in 2nd direction, 6 in z-direction 
    else:
        print('Plese provide correct cell type.')
        return 0
    

    # grid where each entry contains contribution from Li and Nb individually for 1 chain
    displ_grid = np.zeros( np.append( dimensions*num_direction, [2, 3] ) ) 


    # ---------------------------------- Here happens the magic now: ----------------------------------

    # iterate over all atoms
    for site in structure.sites: 
        species = site.specie.symbol

        # only calculate contribution, if atom is Li, Nb, or Ta (we define polarization such that oxygen has no contribution)
        if species in atoms_dict.keys():

            # get assigned charge + assigned number of nearest neighbors
            # at this point, we could also take value of BEC instead of nominal charge....!
            chg, num_nn, idx = atoms_dict[species] 

            # get all neighbors in a sphere of radius nn_dist -> need this to calculate the displacement from paraelectric site
            all_neighbors = structure.get_neighbors( site, nn_dist )
            if len(all_neighbors) < 6:
                    print(f'error with {species}! Only found {len(all_neighbors)} neighbors')
            
            # determine nearest neighbors. 3 for Li, 6 for Nb/Ta 
            nn_mask = np.argsort( [ abs(s2.frac_coords[2] - site.frac_coords[2]) for s2 in all_neighbors] )
            # find center of charge for oxygen plane/cage -> just the mean of nearest neighbors
            center_of_oxygen_charge = np.mean( [all_neighbors[i].frac_coords for i in nn_mask[0:num_nn]], axis=0 )

            # get displacement of cation from praelectric lattice site (i.e. oxygen center of charge) -> this is a 3D vector in angstroems!
            displ = (site.frac_coords - center_of_oxygen_charge) * structure.lattice.abc

            # now we got dipole moment induced by one cation. But where in our grid does it sit?
            # determine coordinates of calcualted dipole moment in the grid of output
            grid_pos = np.array( np.rint(site.frac_coords[0:3] * num_direction * dimensions), dtype=int )
            grid_pos = np.mod( grid_pos, dimensions*num_direction )
            displ_grid[grid_pos[0], grid_pos[1], grid_pos[2], idx] += displ


    # renormalize dipole moments:
    layer_volume = structure[0].lattice.volume / (dimensions[args.normalaxis] * num_direction[args.normalaxis])
    factor = 1e3*1.6/layer_volume # with this factor we get polarization in muC/cm^2
    displ_grid *= factor 

    # list of coordinates normal to wall
    x_n = np.linspace(0, structure[0].lattice.abc[args.normalaxis], num_direction[args.normalaxis]*dimensions[args.normalaxis], endpoint=False)

    # obtain the two axes parallel to the wall
    wall_axes = tuple( np.delete( [0, 1, 2], args.normalaxis ) )
    
    # Calculate averaged polarizations:
    # z-component; 1: Li contribuion; 2: Nb/Ta contribution:
    p_1_z = np.sum(displ_grid[:, :, :, 0, 2], axis=wall_axes) * atoms_dict['Li'][0] # Li contribution
    p_2_z = np.sum(displ_grid[:, :, :, 1, 2], axis=wall_axes) * atoms_dict[args.system][0] # Nb/Ta contribution
    # x-component:
    p_1_x = np.sum(displ_grid[:, :, :, 0, 0], axis=wall_axes) * atoms_dict['Li'][0]
    p_2_x = np.sum(displ_grid[:, :, :, 1, 0], axis=wall_axes) * atoms_dict[args.system][0]
    # y-component:
    p_1_y = np.sum(displ_grid[:, :, :, 0, 1], axis=wall_axes) * atoms_dict['Li'][0]
    p_2_y = np.sum(displ_grid[:, :, :, 1, 1], axis=wall_axes) * atoms_dict[args.system][0]
    
    # ---------------------------------- Save results: ----------------------------------
    outdata = np.array([x_n, p_1_x, p_1_y, p_1_z, p_2_x, p_2_y, p_2_z])
    header = f' {args.xmult}x{args.xmult}x1 supercell.  \nx_n     p_Li,x      p_Li,y     p_Li,z     p_{args.system},x     p_{args.system},y     p_{args.system},z'

    out_file = ''
    if args.xmult > args.ymult:
        out_file = args.outfile + f'_{args.xmult}x.dat'
    elif args.ymult > args.xmult:
        out_file = args.outfile + f'_{args.ymult}y.dat'
    else:
        out_file = args.outfile + '.dat'
    
    np.savetxt(out_file, outdata.T, fmt='%.6e', header=header)

    print(f'--------------------------------------\nExtraced all polarizations; wrote to {out_file}\n--------------------------------------\n')

    # plot the curves if necessacy for debug purposes
    if args.plot_pz > 0:
         import matplotlib.pyplot as plt

         plt.plot(x_n, p_1_z)
         plt.plot(x_n, p_2_z)
         plt.plot(x_n, p_1_z+p_2_z)
         plt.show()



if __name__ == "__main__":
    main()