<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D CG for Hex</description>
    <executable>Helmholtz3D</executable>
    <parameters>helmholtz3D_hex.xml</parameters>
    <files>
        <file description="Session File">helmholtz3D_hex.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.11818e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">5.66374e-07</value>
        </metric>
    </metrics>
</test>


