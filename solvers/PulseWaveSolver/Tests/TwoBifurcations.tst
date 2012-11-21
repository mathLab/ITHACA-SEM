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
            <value variable="A" tolerance="1e-12">43.9048</value>
            <value variable="u" tolerance="1e-12">0.0152826</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">9.81818</value>
            <value variable="u" tolerance="1e-12">0.0100346</value>
        </metric>
    </metrics>
</test>


