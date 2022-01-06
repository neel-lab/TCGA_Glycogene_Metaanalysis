#!/bin/bash

#1. Parse the Ontology into a list datastructure:
Rscript ./glycoOnto_Graph.R
#2. Create Pathway DF:
Rscript ./create_pathLists.R
