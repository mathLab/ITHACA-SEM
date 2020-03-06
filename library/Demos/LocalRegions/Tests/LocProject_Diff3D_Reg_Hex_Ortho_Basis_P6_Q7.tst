<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff3D Reg. Hex Orthogonal Basis, P=6, Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s hexahedron -b Ortho_A Ortho_A Ortho_A -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1.0 0.0 1.0 1.0 1.0 1.0 0.0 1.0 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.72684e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-08">2.36469e-11</value>
        </metric>
    </metrics>
</test>


