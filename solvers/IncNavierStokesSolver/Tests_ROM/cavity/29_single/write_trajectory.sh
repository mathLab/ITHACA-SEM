

#!/bin/sh
for i in {58415..59414}
do
  FieldConvert cavity_poi_29_stage1.xml cavity_poi_29_stage1_$i.chk  mov_$i.vtu
done



