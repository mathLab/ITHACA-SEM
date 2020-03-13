<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>55 Artery Network, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>55_Artery_Network.xml</parameters>
    <files>
        <file description="Session File">55_Artery_Network.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">35.2833</value>
            <value variable="u" tolerance="1e-12">66.3241</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">6.35242</value>
            <value variable="u" tolerance="1e-12">22.2798</value>
        </metric>
    </metrics>
</test>


