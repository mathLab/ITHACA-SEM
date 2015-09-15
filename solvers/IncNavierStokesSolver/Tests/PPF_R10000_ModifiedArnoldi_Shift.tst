<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (Modified Arnoldi with shift): ChannelMax Ev = (3.7302e-03+2.3752e-01i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>PPF_R10000_3D.xml</parameters>
    <files>
        <file description="Session File">PPF_R10000_3D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="2">
            <value variable="u" tolerance="0.0001">3.01846</value>
            <value variable="v" tolerance="0.0001">1.8469</value>
            <value variable="w" tolerance="0.0001">5.31914e-06</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="0.0001">2.25648</value>
            <value variable="v" tolerance="0.0001">0.985804</value>
            <value variable="w" tolerance="0.0001">1.08038e-05</value>
        </metric>
    </metrics>
</test>


