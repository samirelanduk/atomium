/*

This code written by Professor Andrew Martin on 2019-08-08

alias cc='/usr/bin/cc -I$HOME/include -L$HOME/lib -ansi -pedantic -Wall'
cc -o readapdb readapdb.c -lbiop -lgen -lm -lxml2
*/

#include <stdio.h>
#include "bioplib/pdb.h"

int main(int argc, char **argv)
{
   FILE *fp;
   
   if((fp=fopen(argv[1], "r"))!=NULL)
   {
      int natoms;
      PDB *pdb;
      
      if((pdb = blReadPDB(fp, &natoms))!=NULL)
      {

         PDBSTRUCT *pdbs;
         pdbs = blAllocPDBStructure(pdb);

      }
   }
   
   return(0);
}

