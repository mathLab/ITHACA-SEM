<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>55 Artery Network, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>Paper4_55Network.xml</parameters>
    <files>
        <file description="Session File">Paper4_55Network.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">12.286</value>
            <value variable="u" tolerance="1e-12">19.2064</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">6.1916</value>
            <value variable="u" tolerance="1e-12">11.8315</value>
        </metric>
    </metrics>
</test>


