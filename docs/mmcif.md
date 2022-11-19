# Guide to mmCIF

### `entry`

A single row table giving the structure's PDB code.

## Entities

### `entity`

Lists the entities in the file - the unique molecule types. These can be 'polymer', 'carbohydrate', 'non-polymer' or 'water'.

**PDB Records**: `COMPND` mostly, as it is the only way to work out which polymer chains are grouped into entities.

### `entity_name_com`

For each polymer entity, a set of alternative names are provided.

**PDB Records**: `COMPND`

### `entity_poly`

For each polymer entity, the canonical sequence is given, as well as the IDs of the polymers that belong to this entity.

**PDB Records**: `SEQRES`

## Coordinates

### `atom_type`

An alphabetical list of elements that appear in the coordinates.

**PDB Records**: `ATOM`, `HETATM`,

### `atom_site`

A list of atoms with their coordinates.

The `asym_id` of an atom is the molecule it belongs to, and there are `auth` and `label` ones. For polymers, the PDB equivalent is used. For non-polymers, the PDB equivalent is used for `auth` (so, the chain it is associated with) and a new ID is used for `label`. For waters, the same rule applies but all waters associated with a particular chain are treated as though they are a single molecule.

The sequence number is taken from the PDB equivalent for all atoms for the `auth` value, but only for polymers for the `label` value.

**PDB Records**: `ATOM`, `HETATM`, `TER`, `MODEL`, `ENDMDL`

### `atom_site_anisotrop`

Anisotropy associated with atoms.

**PDB Records**: `ANISOU`