<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 3D CG for Pyramid</description>
    <executable>Helmholtz3D</executable>
    <parameters>Helmholtz3D_Pyr.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Pyr.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.00954e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">4.36308e-06</value>
        </metric>
    </metrics>
</test>
