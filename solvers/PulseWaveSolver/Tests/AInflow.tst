<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>AInflow Inlet Condition, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>AInflow.xml</parameters>
    <files>
        <file description="Session File">AInflow.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">14.7991</value>
            <value variable="u" tolerance="1e-12">15.4943</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-4">1.01308</value>
            <value variable="u" tolerance="1e-4">6.50395</value>
        </metric>
    </metrics>
</test>
