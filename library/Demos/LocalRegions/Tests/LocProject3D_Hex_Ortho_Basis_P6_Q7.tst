<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project3D Hex Ortho basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s hexahedron -b Ortho_A Ortho_A Ortho_A -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1.0 0.0 1.0 1.0 1.0 1.0 0.0 1.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.41906e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">5.11591e-13</value>
        </metric>
    </metrics>
</test>
