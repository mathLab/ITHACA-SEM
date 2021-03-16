

#!/bin/sh
for i in {63773..66772}
do
  FieldConvert cavity_poi_1000_stage1.xml cavity_poi_1000_stage1_$i.chk  mov_$i.vtu
done



