<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Double Bifurcation, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>TwoBifurcations.xml</parameters>
    <files>
        <file description="Session File">TwoBifurcations.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">96.2931</value>
            <value variable="u" tolerance="1e-12">11.9121</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">6.31075</value>
            <value variable="u" tolerance="1e-12">6.15717</value>
        </metric>
    </metrics>
</test>


