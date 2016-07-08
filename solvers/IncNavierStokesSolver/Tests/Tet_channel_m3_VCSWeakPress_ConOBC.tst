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
            <value variable="u" tolerance="1e-12">3.41861e-09</value>
            <value variable="v" tolerance="1e-12">6.63202e-09</value>
            <value variable="w" tolerance="1e-12">9.29286e-09</value>
	    <value variable="p" tolerance="1e-12">2.72933e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.27905e-08</value>
            <value variable="v" tolerance="1e-12">3.08842e-08</value>
            <value variable="w" tolerance="1e-12">5.9342e-08</value>
	    <value variable="p" tolerance="1e-12">2.90078e-05</value>
        </metric>
    </metrics>
</test>

