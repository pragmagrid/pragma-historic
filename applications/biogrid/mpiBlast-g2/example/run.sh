#!/bin/sh
MPIBLAST_HOME=$HOME/pkg/mpiblast_g2
$MPIBLAST_HOME/bin/blast-g2.submit -compressio -debug -prof -program blastn -query BLAST.in -database drosoph.nt -machinefile machinefile &
