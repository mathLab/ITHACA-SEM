<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject3D Tetrahedron Modified basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s tetrahedron -b Modified_A Modified_B Modified_C -o 6 6 6 -p 7 7 7</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.89776e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="2e-12">4.32543e-13</value>
        </metric>
    </metrics>
</test>


