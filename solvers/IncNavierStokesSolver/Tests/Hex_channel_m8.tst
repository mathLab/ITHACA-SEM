<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hexahedral elements, P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_channel_m8.xml</parameters>
    <files>
        <file description="Session File">Hex_channel_m8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">7.86396e-14</value>
            <value variable="v" tolerance="1e-8">1.3903e-13</value>
            <value variable="w" tolerance="1e-8">4.50797e-13</value>
	    <value variable="p" tolerance="1e-8">1.18582e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">3.06127e-13</value>
            <value variable="v" tolerance="1e-8">4.79431e-13</value>
            <value variable="w" tolerance="1e-8">1.53633e-12</value>
	    <value variable="p" tolerance="1e-8">2.57923e-11</value>
        </metric>
    </metrics>
</test>
