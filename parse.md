## mmCIF Dictionary to File Object

The entities in the model are determined from the `_entity` table, and annotated with extra information from related tables (such as sequence information, in the case of polymer entities). These are the types of molecule that will be found, with each molecule type being either a polymer, branched chain, non-polymer, or water. Each is turned into a dynamically generated class, which subclasses either `Polymer`, `BranchedPolymer` `NonPolymer`,  or `Water`.

The secondary structure tables (`_struct_conf` and `_struct_sheet_range`) are then checked to get helix and strand information - each of which is a list of lists of residue IDs.

The `_atom_site` table is then split into the different models. Each of the following steps is applied to these individual lists of atoms.

The list of discrete structures in the model are determined by going through the atoms and looking at the IDs there. Each of these structures has a label_asym_id and an auth_asym_id, and is linked to one of the entities. atomium expects that each structure's entity ID will match one in the `_entity` table.

Each of these lists of atoms is turned into a molecule object (specifically an instance of one of the entity classes). These are then passed to a Model object, and the models collectively are passed to the File.