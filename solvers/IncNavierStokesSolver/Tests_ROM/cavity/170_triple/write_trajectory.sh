

#!/bin/sh
for i in {61771..62770}
do
  FieldConvert cavity_poi_170_stage2.xml cavity_poi_170_stage2_$i.chk  mov_$i.vtu
done



