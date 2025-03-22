# Add NME and ACE chemical groups to protein terminal residues
## Description
This script adds NME and ACE groups to protein terminal residues. The pdb file may contain one or more than one chains. In case of multiple chains, the scripts adds NME and ACE to each chains. Note that chains are detected from the chain labels like 'A', 'B' etc using MDAnalysis `mda.segments` function. Also make sure that terminal residues have backbone `C`, `O` and `CA` atoms. Also note that the output **removes TER** in multi chain systems, so before running tleap add TER between chains. It can be done with a quick bash command that adds TER between lines where one line has NME and the next one has ACE, e.g something like

        awk '/NME/{nme=NR} /ACE/ && nme && NR > nme {print "TER"; nme=0} {print}' capped.pdb > capped_TER.pdb

## Dependencies
1. Python3 (I tested with 3.11.7), should also work for other python3 versions
2. MDAnalysis (I tested with 2.7.0), should also work for other versions
3. Numpy (I tested with 1.26.4), but should also work with other versions

Before running the command make sure that the input protein has all the waters and other stuff removed. Also **remove all the Hydrogens**, tleap can add it automatically later. You can remove the hydrogen by using the following quick command, considering your protein is named `protein.pdb`:

    python -c 'import MDAnalysis as mda; u=mda.Universe ("protein.pdb"); u.select_atoms("protein and not name H* 1H* 2H* 3H*").write ("input.pdb")'

## Running the command

Now to add ACE and NME to the hydrogen-removed-pdb, run following command:

         python add_caps.py -i input.pdb -o output.pdb
         
## How it works
  The script reads the hydrogen-less pdb file using MDAnalysis. 
### To add NME
1. If the first residue has an `OXT`, we remove it and place an `N` atom of NME. The carbon atom of NME is added at a distance of 1.36 Å along the `C-N` vector. Where `C` is the protein backbone carbon atom.
2. If there is no `OXT`,  the `N` atom is connected to the backbone `C` along the vector connecting the backbone `C` position and the mid-point of `CA` and `O` atoms of the backbone.
### To add ACE
1. The first carbon of ACE is connected to the backbone `N` and along the vector connecting backbone `CA` and `N` atoms.
2. To add the other carbon and oxygen of ACE, we imagine the backbone `N`, the other carbon and oxygen to be the vertices of an equilateral triangle with the previously added carbon as the centroid of the triangle. The coordinates of other two vertices given one vertex ($x_a, y_a, z_a$) and centroid ($x_g, y_g, z_g$) of an equilateratl triangle is obtained as:

   
  First vertex:
  
  
  $x_1 = x_g - \frac{(x_a - x_g)}{2} + \frac{\sqrt{3}}{2} \left( n_y (z_a - z_g) - n_z (y_a - y_g) \right)$
  
  
  $y_1 = y_g - \frac{(y_a - y_g)}{2} + \frac{\sqrt{3}}{2} \left( n_z (x_a - x_g) - n_x (z_a - z_g) \right)$
  
  $z_1 = z_g - \frac{(z_a - z_g)}{2} + \frac{\sqrt{3}}{2} \left( n_x (y_a - y_g) - n_y (x_a - x_g) \right)$
  
  Second vertex:
  
  $x_2 = x_g - \frac{(x_a - x_g)}{2} - \frac{\sqrt{3}}{2} \left( n_y (z_a - z_g) - n_z (y_a - y_g) \right)$
  
  $y_2 = y_g - \frac{(y_a - y_g)}{2} - \frac{\sqrt{3}}{2} \left( n_z (x_a - x_g) - n_x (z_a - z_g) \right)$
  
  $z_2 = z_g - \frac{(z_a - z_g)}{2} - \frac{\sqrt{3}}{2} \left( n_x (y_a - y_g) - n_y (x_a - x_g) \right)$
  
  
  Where ($n_x, n_y, n_z$) is an arbritrary orientation as in 3D the two vertices can be rotated about the other vertex and the centroid abritrarily i.e infinite solutions. We select these unit vector randomly using `np.random`, it really does not matter. The coordinates generated above may be a over 2 Å apart. So finally we rescale the bonds such that they are around 1.4 Å. 

### Create universe
MDAnalysis universe are created for above `ACE` and `NME`. 
  1. The names of `ACE` atoms are **C, CH3, O**. While the NME atoms are named **N, C**. Note that in earlier ambertools the NME are named **N, CH3**, so one can replace C with CH3 for such cases if later on tleap complains about atom names not being found.
 2. The ACE, protein and NME universes are merged with MDAnalysis `mda.Merge` function and the final pdb is written
