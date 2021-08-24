<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project3D Deformed Hex Lagrange basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s hexahedron -b GLL_Lagrange GLL_Lagrange GLL_Lagrange -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.5 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1.5 0.0 1.0 1.0 1.0 1.0 0.0 1.0 1.5</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.00820104</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">0.0804162</value>
        </metric>
    </metrics>
</test>
