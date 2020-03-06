<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff3D Reg. Tet Modified Basis, P=6, Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s tetrahedron -b Modified_A Modified_B Modified_C -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">9.02964e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">5.08749e-12</value>
        </metric>
    </metrics>
</test>


