

#!/bin/sh
for i in {0..1000}
do
  FieldConvert ICOSAHOM2018_compare_RB_sparsePoly.xml ICOSAHOM2018_compare_RB_sparsePoly_$i.chk  mov_$i.vtu
done



