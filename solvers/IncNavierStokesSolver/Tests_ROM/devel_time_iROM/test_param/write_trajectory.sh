

#!/bin/sh
for i in {0..25}
do
  FieldConvert ICOSAHOM2018_compare_RB_sparsePoly.xml ICOSAHOM2018_compare_RB_sparsePoly_$i.chk  mov_$i.vtu
done

for i in {0..25}
do
  FieldConvert ICOSAHOM2018_compare_RB_sparsePoly_ROM2.xml ICOSAHOM2018_compare_RB_sparsePoly_ROM2_$i.chk  movROM_$i.vtu
done


