<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (Modified Arnoldi with shift): ChannelMax Ev = (3.7302e-03+2.3752e-01i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>PPF_R10000_3D.xml</parameters>
    <files>
        <file description="Session File">PPF_R10000_3D.xml</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value tolerance="0.001">0.518448,-26.6405</value>
            <value tolerance="0.001">0.518448,26.6405</value>
        </metric>
    </metrics>
</test>


