<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project3D Hex Modified basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s hexahedron -b Modified_A Modified_A Modified_A -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1.0 0.0 1.0 1.0 1.0 1.0 0.0 1.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-08">5.65299e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-08">1.74836e-10</value>
        </metric>
    </metrics>
</test>
