<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff3D Lin. Deformed Hex Lagrange Basis, P=6, Q=7</description>
    <executable>LocProject_Diff3D</executable>
    <parameters>8 8 8 8 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.0353098</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0.446618</value>
        </metric>
    </metrics>
</test>


