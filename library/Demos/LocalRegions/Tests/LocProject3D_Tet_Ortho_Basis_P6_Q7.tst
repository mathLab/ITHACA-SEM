<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project3D Tet Ortho basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s tetrahedron -b Ortho_A Ortho_B Ortho_C -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">4.91242e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">7.99361e-15</value>
        </metric>
    </metrics>
</test>
