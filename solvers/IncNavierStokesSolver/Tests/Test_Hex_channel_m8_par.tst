<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hex elements, par(2), P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Test_Hex_channel_m8_par.xml</parameters>
    <files>
        <file description="Session File">Test_Hex_channel_m8_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.03412e-11</value>
            <value variable="v" tolerance="1e-12">7.62679e-12</value>
            <value variable="w" tolerance="1e-12">8.14841e-11</value>
	    <value variable="p" tolerance="1e-12">2.21238e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">9.31369e-11</value>
            <value variable="v" tolerance="1e-12">7.31505e-11</value>
            <value variable="w" tolerance="1e-12">5.83095e-09</value>
	    <value variable="p" tolerance="1e-12">9.42423e-08</value>
        </metric>
    </metrics>
</test>