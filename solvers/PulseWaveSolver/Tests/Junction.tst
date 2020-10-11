<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Junction Riemann Solver, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>Junction.xml</parameters>
    <files>
        <file description="Session File">Junction.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">79.0903</value>
            <value variable="u" tolerance="1e-12">5.63277e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">5.983</value>
            <value variable="u" tolerance="1e-12">4.25274e-06</value>
        </metric>
    </metrics>
</test>


