

#!/bin/sh
for i in {61770..62769}
do
  FieldConvert cavity_poi_150_stage2.xml cavity_poi_150_stage2_$i.chk  mov_$i.vtu
done



