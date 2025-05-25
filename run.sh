#!/bin/bash

for N in 1000 2000 3000 4000
do
  echo "Generating input for N=$N"
  ./generate_input $N input${N}.txt
done
