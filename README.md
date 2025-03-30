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
The NME `N` atom is added to satisfy a trigonal planar geometry around terminal residue `C` using the `C`, `O`, and `CA` positions as inputs. The NME `C`, is added to satisfy the NME `N`-`C` bond length, the residue `C`, NME `N`, NME `C` angle of 120deg, and residue `O`, residue `C`, NME `N`, NME `C` dihedral of 0deg.

### To add ACE
ACE `C`, `CH3`, and `O` atoms are added using the expected bond length, bond angles, and dihedral angles implied by the backbone of the N-terminal residue.


### Create universe
MDAnalysis universe are created for above `ACE` and `NME`. 
  1. The names of `ACE` atoms are **C, CH3, O**. While the NME atoms are named **N, C**. Note that in earlier ambertools the NME are named **N, CH3**, so one can replace C with CH3 for such cases if later on tleap complains about atom names not being found.
 2. The ACE, protein and NME universes are merged with MDAnalysis `mda.Merge` function and the final pdb is written
