<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (LM with Arpack): ChannelMax Ev = (0.0037303,+/-0.237523i) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>PPF_R10000.xml</parameters>
    <files>
        <file description="Session File">PPF_R10000.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="2">
            <value variable="u" tolerance="0.00001">0.0409899</value>
            <value variable="v" tolerance="0.00001">0.0168088</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="0.0001">0.0275149</value>
            <value variable="v" tolerance="0.0001">0.0119663</value>
        </metric>
    </metrics>
</test>


