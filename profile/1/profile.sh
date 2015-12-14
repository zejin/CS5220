#!/bin/sh

module load maqao

amplxe-cl -report hotspots -r r000ah/ > amplxe_all.txt
amplxe-cl -report hotspots -source-object function="main" > amplxe_main.txt
maqao cqa ./sim fct=main uarch=HASWELL > maqao.txt