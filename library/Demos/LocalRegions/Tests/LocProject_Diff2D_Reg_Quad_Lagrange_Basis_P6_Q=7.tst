<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff2D Reg. Quad Lagrange Basis P=6, Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s quadrilateral -b GLL_Lagrange GLL_Lagrange -o 6 6 -p 7 7 -c 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">6.80298e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">4.83169e-13</value>
        </metric>
    </metrics>
</test>


