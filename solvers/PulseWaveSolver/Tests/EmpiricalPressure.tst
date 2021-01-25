<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Empirical PressureArea Law, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>EmpiricalPressure.xml</parameters>
    <files>
        <file description="Session File">EmpiricalPressure.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">14.807</value>
            <value variable="u" tolerance="1e-12">15.1812</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">1.01605</value>
            <value variable="u" tolerance="1e-12">4.87762</value>
        </metric>
    </metrics>
</test>
