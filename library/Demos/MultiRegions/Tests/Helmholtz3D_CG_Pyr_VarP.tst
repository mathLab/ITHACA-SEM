<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 3D CG for Pyramid with Variable P</description>
    <executable>Helmholtz3D</executable>
    <parameters>Helmholtz3D_Pyr_VarP.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Pyr_VarP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-10">7.06911e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">0.000510411</value>
        </metric>
    </metrics>
</test>
