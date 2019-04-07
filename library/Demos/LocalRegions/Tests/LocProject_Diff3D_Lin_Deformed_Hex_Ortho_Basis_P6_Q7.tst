<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff3D Lin. Deformed Hex Ortho Basis, P=6, Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s hexahedron -b Ortho_A Ortho_A Ortho_A -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.5 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1.5 0.0 1.0 1.0 1.0 1.0 0.0 1.0 1.5 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.0353098</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0.446618</value>
        </metric>
    </metrics>
</test>


