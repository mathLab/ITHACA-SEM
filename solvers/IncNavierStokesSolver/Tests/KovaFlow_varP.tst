<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow variable P</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_varP.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_varP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">8.95241e-05</value>
            <value variable="v" tolerance="1e-12">6.41405e-06</value>
	    <value variable="p" tolerance="1e-12">0.000932063</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000144453</value>
            <value variable="v" tolerance="1e-12">9.04985e-06</value>
	    <value variable="p" tolerance="1e-12">0.00283058</value>
        </metric>
    </metrics>
</test>
