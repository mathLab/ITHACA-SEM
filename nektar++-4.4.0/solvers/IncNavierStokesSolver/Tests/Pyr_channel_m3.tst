<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Pyramidic elements, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Pyr_channel_m3.xml</parameters>
    <files>
        <file description="Session File">Pyr_channel_m3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-8">1.16849e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.75668e-16</value>
            <value variable="v" tolerance="1e-12">1.07627e-16</value>
            <value variable="w" tolerance="1e-12">5.55112e-16</value>
            <value variable="p" tolerance="1e-8">4.04121e-14</value>
        </metric>
    </metrics>
</test>
