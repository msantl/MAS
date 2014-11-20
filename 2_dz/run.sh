#! /bin/bash

for ((i=0; i < 1024; ++i))
do 
    echo -ne "$i    "; ./matija_santl_dz2.exe $i
done > out.txt
