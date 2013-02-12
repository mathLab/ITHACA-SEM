<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hexahedral elements, P=8, Successive RHS(5)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_channel_m8_srhs.xml</parameters>
    <files>
        <file description="Session File">Hex_channel_m8_srhs.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">2.73101e-11</value>
            <value variable="v" tolerance="1e-8">2.23013e-11</value>
            <value variable="w" tolerance="1e-8">6.29721e-10</value>
            <value variable="p" tolerance="1e-8">1.94395e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">1.55008e-10</value>
            <value variable="v" tolerance="1e-8">1.58963e-10</value>
            <value variable="w" tolerance="1e-8">8.27366e-08</value>
            <value variable="p" tolerance="1e-8">2.42343e-08</value>
        </metric>
    </metrics>
</test>
