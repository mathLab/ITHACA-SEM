

#!/bin/sh
for i in {60770..61769}
do
  FieldConvert cavity_poi_150_stage1.xml cavity_poi_150_stage1_$i.chk  mov_$i.vtu
done



