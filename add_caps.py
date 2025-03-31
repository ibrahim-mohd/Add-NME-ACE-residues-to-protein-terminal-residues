# Written by Mohd Ibrahim
# Technical University of Munich
# Email: ibrahim.mohd@tum.de
import numpy as np
import MDAnalysis as mda
import argparse
import warnings
# Suppress specific warnings from MDAnalysis
warnings.filterwarnings("ignore")#, category=UserWarning, module="MDAnalysis.coordinates.PDB")


parser = argparse.ArgumentParser(description="Add capping groups ACE and NME to protein termini. Remove the hydrogens from the input pdb file before using this script")
parser.add_argument('-i', dest='in_file', type=str, default='protein_noh.pdb',help='pdb file')
parser.add_argument('-o', dest='out_file', type=str, default='protein_noh_cap.pdb',help='output file')

args      = parser.parse_args()
in_file   = args.in_file
out_file  = args.out_file

def create_universe (n_atoms, name, resname, positions, resids, segid):

    u_new = mda.Universe.empty(n_atoms=n_atoms,
                             n_residues=n_atoms,
                             atom_resindex=np.arange (n_atoms),
                             residue_segindex=np.arange (n_atoms),
                             n_segments=n_atoms,
                             trajectory=True) # necessary for adding coordinate


    u_new.add_TopologyAttr('name',   name)
    u_new.add_TopologyAttr('resid', resids)
    u_new.add_TopologyAttr('resname', resname)
    u_new.atoms.positions = positions
    u_new.add_TopologyAttr('segid', n_atoms*[segid])
    u_new.add_TopologyAttr('chainID', n_atoms * [segid])


    return u_new

def calcCoordinate(a, b, c, bond_len, theta, di_angle):
    """
    Calculate position of fourth atom 'd' based on input atoms a,b,c
    and known geometry wrt a,b,c.
    bond_len: c-d bond length (angstrom)
    theta: b,c,d angle (degrees)
    di_angle: a,b,c,d dihedral (degrees)

    with thanks to https://github.com/osita-sunday-nnyigide/Pras_Server/
    """
    di_angle = np.deg2rad(di_angle)

    u = c - b
    x = a - b
    v = ((x) - (np.dot((x), u)/np.dot(u,u))*u)

    w = np.cross(u, x)

    q = (v/np.linalg.norm(v))*np.cos(di_angle)
    e = (w/np.linalg.norm(w))*np.sin(di_angle)

    pos_temp2 = np.array((b + (q+e)))

    u1 = b - c
    y1 = pos_temp2 - c

    mag_y1 = np.linalg.norm(y1)
    mag_u1 = np.linalg.norm(u1)

    theta_bcd = np.arccos(np.dot(u1,y1)/(mag_u1*mag_y1))
    rotate = np.deg2rad(theta) - theta_bcd

    z  = np.cross(u1, y1)
    n = z/np.linalg.norm(z)

    pos_ini = c + y1*np.cos(rotate) +\
              (np.cross(n,y1)*np.sin(rotate))\
              + n*(np.dot(n,y1))*(1-np.cos(rotate))

    d = ((pos_ini-c)*(bond_len/np.linalg.norm(pos_ini-c)))+c

    return d

def get_nme_pos (end_residue):

    # find midpoint of O and CA
    index_o = np.where (end_residue.names == "O")[0][0]
    index_ca = np.where (end_residue.names == "CA")[0][0]
    index_c = np.where(end_residue.names=='C')[0][0]

    pos_o = end_residue.positions[index_o]
    pos_ca = end_residue.positions[index_ca]
    pos_c = end_residue.positions[index_c]

    #bisect the O, C, CA angle
    v1 = pos_o-pos_c
    v1 /= np.linalg.norm(v1)
    v2 = pos_ca-pos_c
    v2 /= np.linalg.norm(v2)
    bisector = v1 + v2
    bisector /= np.linalg.norm(bisector)
    # apply translation to NME N atom from main chain C atom
    # in direction opposite bisector:
    bondLength = 1.34
    N_position = bondLength * -(bisector) + pos_c

    # now place NME [C] position using bond length, angle,
    # and dihedral angle of 0deg, where dihedral is O=C-N-[C]
    bondLength = 1.45
    C_position = calcCoordinate(pos_o, pos_c, N_position, bondLength, 120, 0)

    return N_position, C_position

def get_ace_pos (end_residue):

    index_c = np.where (end_residue.names == "C")[0][0]
    index_ca = np.where (end_residue.names == "CA")[0][0]
    index_n  = np.where (end_residue.names == "N")[0][0]

    pos_c = end_residue.positions[index_c]
    pos_ca = end_residue.positions[index_ca]
    pos_n = end_residue.positions[index_n]

    # C-CA-N-C dihedral has two minima at -60,60
    C_position = calcCoordinate(pos_c, pos_ca, pos_n, 1.34, 120, -60)
    CH3_position = calcCoordinate(pos_ca, pos_n, C_position, 1.52, 120, 180)
    O_position = calcCoordinate(pos_ca, pos_n, C_position, 1.23, 120, 0)

    return C_position, CH3_position, O_position


##
# Load pdb file
u = mda.Universe (in_file)

# Access each fragment separately
res_start = 0
segment_universes = []

for seg in u.segments:

    chain = u.select_atoms(f"segid {seg.segid}")

    # Add ACE
    resid_c = chain.residues.resids [0]
    end_residue = u.select_atoms(f"segid {seg.segid} and resid {resid_c}")
    ace_positions = get_ace_pos (end_residue)
    ace_names = ["C", "CH3", "O"]
    resid = chain.residues.resids[0]
    kwargs = dict (n_atoms=len(ace_positions), name=ace_names,
                      resname=len(ace_names)*["ACE"], positions=ace_positions,
                      resids=resid*np.ones(len(ace_names)),
                      segid=chain.segids[0])

    ace_universe =  create_universe (**kwargs)

    # Add NME
    resid_c     = chain.residues.resids [-1]
    end_residue = u.select_atoms(f"segid {seg.segid} and resid {resid_c}")


    nme_positions = get_nme_pos (end_residue)
    nme_names   = ["N", "C"]

    resid = chain.residues.resids[-1]+2

    kwargs = dict (n_atoms=len(nme_names), name=nme_names,
                      resname=len(nme_names)*["NME"], positions=nme_positions,
                      resids=resid*np.ones(len(nme_names)),
                      segid=chain.segids[0])

    nme_universe =  create_universe (**kwargs)
    ## Merge Universe
    if "OXT" in end_residue.names:

        index = np.where (end_residue.names == "OXT")[0][0]
        OXT   = end_residue [index]

        Chain     = u.select_atoms(f"segid {seg.segid} and not index {OXT.index}")

    else:

        Chain     = u.select_atoms(f"segid {seg.segid}")

    ### Merge ACE, Protien and NME

    u_all = mda.Merge (ace_universe.atoms, Chain, nme_universe.atoms)


    # to renumber residues
    resids_ace = [res_start+1, res_start+1, res_start+1]
    resids_pro = np.arange (resids_ace[0]+1, Chain.residues.n_residues+resids_ace[0]+1)
    resids_nme = [resids_pro[-1]+1,resids_pro[-1]+1]

    u_all.atoms.residues.resids =  np.concatenate ([resids_ace,resids_pro,resids_nme])#np.arange (1+res_start, len(u_all.atoms.residues.resids)+res_start+1)

    res_start = u_all.atoms.residues.resids[-1]

    segment_universes.append (u_all)

## Join all the universes
all_uni = mda.Merge(*(seg.atoms for seg in segment_universes))

all_uni.atoms.write (out_file)
