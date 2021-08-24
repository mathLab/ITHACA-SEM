<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject_Diff3D Tetrahedron Modified basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s tetrahedron -b Modified_A Modified_B Modified_C -o 6 6 6 -p 7 7 7 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-11">2.1782e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">5.08749e-12</value>
        </metric>
    </metrics>
</test>


