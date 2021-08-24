<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project3D Tet Modified basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s tetrahedron -b Modified_A Modified_B Modified_C -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.89448e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">1.17817e-12</value>
        </metric>
    </metrics>
</test>
