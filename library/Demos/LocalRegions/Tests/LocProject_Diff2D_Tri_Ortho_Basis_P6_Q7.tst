<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff2D Tet Ortho Basis, P=6, Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s triangle -b Ortho_A Ortho_B -o 6 6 -p 7 7 -c 0.0 0.0 1.0 1.0 0.5 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">3.04703e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">5.82645e-13</value>
        </metric>
    </metrics>
</test>


