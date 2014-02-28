<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D pulse with base flow, P=5</description>
    <executable>APESolver</executable>
    <parameters>Pulse.xml</parameters>
    <files>
        <file description="Session File">Pulse.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-12">21.9512</value>
            <value variable="u" tolerance="1e-12">0.00517265</value>
            <value variable="v" tolerance="1e-12">0.00517265</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-12">98.1507</value>
            <value variable="u" tolerance="1e-12">0.0194736</value>
            <value variable="v" tolerance="1e-12">0.0194736</value>
        </metric>
    </metrics>
</test>
