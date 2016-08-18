<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Hexahedral elements, P=8, flowrate</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_channel_m8_Flowrate.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Hex_channel_m8_Flowrate.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.99699e-10</value>
            <value variable="v" tolerance="1e-12">1.56686e-11</value>
            <value variable="w" tolerance="1e-12">1.87622e-09</value>
            <value variable="p" tolerance="1e-8">5.25319e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.04108e-10</value>
            <value variable="v" tolerance="1e-12">2.89958e-11</value>
            <value variable="w" tolerance="1e-12">6.66162e-09</value>
            <value variable="p" tolerance="1e-8">8.27375e-08</value>
        </metric>
    </metrics>
</test>
