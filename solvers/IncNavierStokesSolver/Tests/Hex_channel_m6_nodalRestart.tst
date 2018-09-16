<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Hexahedral elements, P=6, restart from nodal field</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_channel_m6_nodalRestart.xml</parameters>
    <files>
        <file description="Session file">Hex_channel_m6_nodalRestart.xml</file>
        <file description="Restart file">Hex_channel_m6_nodalRestart.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.28479e-14</value>
            <value variable="v" tolerance="1e-12">3.07108e-14</value>
            <value variable="w" tolerance="1e-12">1.17705e-13</value>
            <value variable="p" tolerance="1e-8">1.25369e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.32143e-13</value>
            <value variable="v" tolerance="1e-12">1.38338e-13</value>
            <value variable="w" tolerance="2e-12">4.14974e-13</value>
            <value variable="p" tolerance="1e-8">3.82983e-12</value>
        </metric>
    </metrics>
</test>
