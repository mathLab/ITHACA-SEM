<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>PowerPressureArea, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>PowerPressure.xml</parameters>
    <files>
        <file description="Session File">PowerPressure.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-6">7.5618</value>
            <value variable="u" tolerance="1e-6">54.6041</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-6">1.76251</value>
            <value variable="u" tolerance="1e-6">23.4226</value>
        </metric>
    </metrics>
</test>
