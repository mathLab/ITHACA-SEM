<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project3D Hex Lagrange basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s hexahedron -b GLL_Lagrange GLL_Lagrange GLL_Lagrange -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1.0 0.0 1.0 1.0 1.0 1.0 0.0 1.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-11">2.04616e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">8.52651e-13</value>
        </metric>
    </metrics>
</test>
