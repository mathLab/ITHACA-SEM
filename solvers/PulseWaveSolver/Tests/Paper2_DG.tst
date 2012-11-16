<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D straight vessel, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>Paper2_DG.xml</parameters>
    <files>
        <file description="Session File">Paper2_DG.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">14.1394</value>
            <value variable="u" tolerance="1e-12">14.1161</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">1.00063</value>
            <value variable="u" tolerance="1e-12">1.00633</value>
        </metric>
    </metrics>
</test>


