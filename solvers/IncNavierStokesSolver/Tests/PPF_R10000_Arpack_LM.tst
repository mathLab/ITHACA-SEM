<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (LM with Arpack): ChannelMax Ev = (0.0037303,+/-0.237523i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>PPF_R10000.xml</parameters>
    <files>
        <file description="Session File">PPF_R10000.xml</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value index="0" tolerance="0.001">-0.00024674,0</value>
            <value index="1" tolerance="0.001">-0.00098696,0</value>
        </metric>
    </metrics>
</test>


