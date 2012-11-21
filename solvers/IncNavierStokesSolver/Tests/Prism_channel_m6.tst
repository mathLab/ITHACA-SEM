<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Prismatic elements, P=6</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Prism_channel_m6.xml</parameters>
    <files>
        <file description="Session File">Prism_channel_m6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">4.44507e-16</value>
	    <value variable="p" tolerance="1e-12">1.52206e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">5.28658e-16</value>
            <value variable="v" tolerance="1e-12">3.37939e-16</value>
            <value variable="w" tolerance="1e-12">1.39888e-14</value>
	    <value variable="p" tolerance="1e-12">9.61453e-14</value>
        </metric>
    </metrics>
</test>
