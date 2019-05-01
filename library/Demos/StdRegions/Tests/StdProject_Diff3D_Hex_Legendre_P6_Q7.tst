<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject_Diff3D Hexahedron Legendre basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s hexahedron -b legendre legendre legendre -o 6 6 6 -p 7 7 7 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">4.88883e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.13687e-11</value>
        </metric>
    </metrics>
</test>


