<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Merging Riemann Solver, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>Merging.xml</parameters>
    <files>
        <file description="Session File">Merging.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">80.117</value>
            <value variable="u" tolerance="1e-12">0.00160204</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">5.98311</value>
            <value variable="u" tolerance="1e-12">0.000665193</value>
        </metric>
    </metrics>
</test>


