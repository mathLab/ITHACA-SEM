<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (LM with Arpack and Complex Shift): ChannelMax Ev = (0.00248682 -0.158348i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>PPF_R15000_3D.xml</parameters>
    <files>
        <file description="Session File">PPF_R15000_3D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="2">
            <value variable="u" tolerance="0.0001">0.836151</value>
            <value variable="v" tolerance="0.0001">0.283706</value>
            <value variable="w" tolerance="0.0001">0.473301</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="0.0001">0.1450440</value>
            <value variable="v" tolerance="0.0001">0.0463947</value>
            <value variable="w" tolerance="0.0001">0.10873</value>
        </metric>
    </metrics>
</test>


