<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff2D Reg. Quad Ortho Basis, P=4, Q=5</description>
    <executable>LocProject</executable>
    <parameters>-s quadrilateral -b Ortho_A Ortho_A -o 6 6 -p 7 7 -c 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">6.43909e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">5.68434e-13</value>
        </metric>
    </metrics>
</test>


