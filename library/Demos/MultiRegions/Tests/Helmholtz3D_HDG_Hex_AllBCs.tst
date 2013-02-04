<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 3D HDG, hexes, mixed BCs</description>
    <executable>HDGHelmholtz3D</executable>
    <parameters>Helmholtz3D_Hex_AllBCs_P6.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Hex_AllBCs_P6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.000417674</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0.00168661</value>
        </metric>
    </metrics>
</test>
