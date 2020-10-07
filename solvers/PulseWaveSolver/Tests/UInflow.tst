<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>UInflow Inlet Condition, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>UInflow.xml</parameters>
    <files>
        <file description="Session File">UInflow.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">14.8323</value>
            <value variable="u" tolerance="1e-12">14.8314</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">1.00008</value>
            <value variable="u" tolerance="1e-12">1.0008</value>
        </metric>
    </metrics>
</test>
