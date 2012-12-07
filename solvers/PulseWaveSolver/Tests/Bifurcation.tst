<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Bifurcation Riemann Solver, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>Bifurcation.xml</parameters>
    <files>
        <file description="Session File">Bifurcation.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">59.8338</value>
            <value variable="u" tolerance="1e-12">0.063162</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">6.00991</value>
            <value variable="u" tolerance="1e-12">0.0691512</value>
        </metric>
    </metrics>
</test>


