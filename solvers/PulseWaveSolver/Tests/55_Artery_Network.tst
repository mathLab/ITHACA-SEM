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
            <value variable="A" tolerance="1e-12">34.7647</value>
            <value variable="u" tolerance="1e-12">80.5877</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">5.98239</value>
            <value variable="u" tolerance="1e-12">8.48982</value>
        </metric>
    </metrics>
</test>


