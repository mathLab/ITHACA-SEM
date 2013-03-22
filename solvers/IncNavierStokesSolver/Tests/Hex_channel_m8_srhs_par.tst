<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hexahedral elements, P=8, Successive RHS(5), par(2)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_channel_m8_srhs.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Hex_channel_m8_srhs.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">4.54205e-13</value>
            <value variable="v" tolerance="1e-8">5.92873e-13</value>
            <value variable="w" tolerance="1e-8">2.56712e-12</value>
            <value variable="p" tolerance="1e-8">8.8739e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">1.89675e-12</value>
            <value variable="v" tolerance="1e-8">2.84934e-12</value>
            <value variable="w" tolerance="1e-8">1.51304e-11</value>
            <value variable="p" tolerance="1e-8">3.4398e-10</value>
        </metric>
    </metrics>
</test>
