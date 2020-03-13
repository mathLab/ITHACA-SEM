<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Tet Prism Channel with Radiation outflow </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_prism_channel_rad.xml</parameters>
    <files>
        <file description="Session File">Tet_prism_channel_rad.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
	    <value variable="w" tolerance="1e-12">1.07969e-14</value>
	    <value variable="p" tolerance="1e-12">1.87067e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.21486e-15</value>
            <value variable="v" tolerance="1e-12">1.49877e-15</value>
	    <value variable="w" tolerance="1e-12">1.07969e-14</value>
	    <value variable="p" tolerance="1e-12">1.43219e-13</value>
        </metric>
    </metrics>
</test>


