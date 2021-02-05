

#!/bin/sh
for i in {59769..60768}
do
  FieldConvert cavity_poi_200_stage1.xml cavity_poi_200_stage1_$i.chk  mov_$i.vtu
done



