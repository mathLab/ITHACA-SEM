<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Tetrahedral elements, P=3 with Convective link outflow BC</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_channel_m3_VCSWeakPress_ConOBC.xml</parameters>
    <files>
        <file description="Session File">Tet_channel_m3_VCSWeakPress_ConOBC.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">1.63577e-08</value>
            <value variable="v" tolerance="1e-7">1.09437e-07</value>
            <value variable="w" tolerance="2e-7">2.12151e-07</value>
	    <value variable="p" tolerance="1e-5">1.85022e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">7.63186e-08</value>
            <value variable="v" tolerance="3.5e-7">5.14798e-07</value>
            <value variable="w" tolerance="2e-6">1.69006e-06</value>
	    <value variable="p" tolerance="1e-5">3.47498e-05</value>
        </metric>
    </metrics>
</test>

