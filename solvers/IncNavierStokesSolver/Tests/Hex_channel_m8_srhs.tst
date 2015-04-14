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
            <value variable="u" tolerance="1e-8">5.45235e-13</value>
            <value variable="v" tolerance="1e-8">5.34748e-13</value>
            <value variable="w" tolerance="1e-8">2.40482e-12</value>
            <value variable="p" tolerance="1e-8">4.71821e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">2.01094e-12</value>
            <value variable="v" tolerance="1e-8">2.31639e-12</value>
            <value variable="w" tolerance="1e-8">1.51304e-11</value>
            <value variable="p" tolerance="1e-8">1.04359e-09</value>
        </metric>
    </metrics>
</test>
