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