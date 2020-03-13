<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>RCR outflow boundary condition</description>
    <executable>PulseWaveSolver</executable>
    <parameters>RCR_boundary.xml</parameters>
    <files>
        <file description="Session File">RCR_boundary.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">5.46359e-05</value>
            <value variable="u" tolerance="1e-12">0.427285</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">8.75268e-05</value>
            <value variable="u" tolerance="1e-12">0.781413</value>
        </metric>
    </metrics>
</test>
