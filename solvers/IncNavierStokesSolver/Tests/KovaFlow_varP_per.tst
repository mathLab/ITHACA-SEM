<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow variable P, periodic BC</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_varP_per.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_varP_per.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">7.24234e-05</value>
            <value variable="v" tolerance="1e-12">0.000118842</value>
	    <value variable="p" tolerance="1e-12">0.000545425</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000136274</value>
            <value variable="v" tolerance="1e-12">0.000329206</value>
	    <value variable="p" tolerance="1e-12">0.00303542</value>
        </metric>
    </metrics>
</test>

