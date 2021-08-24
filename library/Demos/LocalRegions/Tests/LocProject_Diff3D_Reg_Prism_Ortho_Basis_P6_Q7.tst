<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff3D Reg. Prism Ortho Basis, P=6, Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s prism -b Ortho_A Ortho_A Ortho_B -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.5 0.0 1.0 0.5 1.0 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-11">2.97892e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">3.12639e-12</value>
        </metric>
    </metrics>
</test>
