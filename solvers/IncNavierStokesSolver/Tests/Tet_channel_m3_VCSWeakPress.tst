<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Tetrahedral elements, Weak Pressure VCS, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_channel_m3_VCSWeakPress.xml</parameters>
    <files>
        <file description="Session File">Tet_channel_m3_VCSWeakPress.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
	    <value variable="p" tolerance="1e-12">1.03457e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.20119e-16</value>
            <value variable="v" tolerance="1e-12">3.48315e-16</value>
            <value variable="w" tolerance="1e-12">.33227e-15</value>
	    <value variable="p" tolerance="1e-12">5.01821e-14</value>
        </metric>
    </metrics>
</test>

