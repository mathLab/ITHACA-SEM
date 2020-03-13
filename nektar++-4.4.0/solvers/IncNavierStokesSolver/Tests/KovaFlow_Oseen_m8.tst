<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_Oseen_m8.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_Oseen_m8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"></value>
            <value variable="v" tolerance="1e-12"></value>
	    <value variable="p" tolerance="1e-12"></value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"></value>
            <value variable="v" tolerance="1e-12"></value>
	    <value variable="p" tolerance="1e-12"></value>
        </metric>
    </metrics>
</test>
