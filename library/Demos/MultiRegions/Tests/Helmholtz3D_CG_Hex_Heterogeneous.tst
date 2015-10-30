<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 3D CG for Hexahedron, heterogeneous expansion</description>
    <executable>Helmholtz3D</executable>
    <parameters>Helmholtz3D_Hex_Heterogeneous.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Hex_Heterogeneous.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.23385e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0.00014837</value>
        </metric>
    </metrics>
</test>
