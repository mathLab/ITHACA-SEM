<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>QInflow Inlet Condition, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>QInflow.xml</parameters>
    <files>
        <file description="Session File">QInflow.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">14.8393</value>
            <value variable="u" tolerance="1e-12">14.9111</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-4">1.05484</value>
            <value variable="u" tolerance="1e-4">1.53557</value>
        </metric>
    </metrics>
</test>
