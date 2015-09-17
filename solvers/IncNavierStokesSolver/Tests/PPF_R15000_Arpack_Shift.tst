<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (LM with Arpack and Complex Shift): ChannelMax Ev = (0.00248682 -0.158348i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>PPF_R15000_3D.xml</parameters>
    <files>
        <file description="Session File">PPF_R15000_3D.xml</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value index="0" tolerance="0.001">-0.0339756,-0.184554</value>
            <value index="1" tolerance="0.001">-0.0339756,-0.215446</value>
        </metric>
    </metrics>
</test>


