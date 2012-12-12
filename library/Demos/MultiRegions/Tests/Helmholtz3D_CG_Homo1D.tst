<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>CG Helmholtz3D Homogeneous 1D</description>
    <executable>Helmholtz3DHomo1D</executable>
    <parameters>Helmholtz3D_Homo1D.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Homo1D.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">0.00423087</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.00918763</value>
        </metric>
    </metrics>
</test>


