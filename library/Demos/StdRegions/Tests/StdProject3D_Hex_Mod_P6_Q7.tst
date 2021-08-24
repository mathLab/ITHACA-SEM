<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject3D Hexahedron Modified basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s hexahedron -b Modified_A Modified_A Modified_A -o 6 6 6 -p 7 7 7</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">4.76767e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">8.0945e-11</value>
        </metric>
    </metrics>
</test>


