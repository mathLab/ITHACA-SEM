<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D CG, hexes, mixed BCs</description>
    <executable>Helmholtz3D</executable>
    <parameters>Helmholtz3D_Hex_MixedBCs_CG.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Hex_MixedBCs_CG.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.000418001</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0.000874515</value>
        </metric>
    </metrics>
</test>


