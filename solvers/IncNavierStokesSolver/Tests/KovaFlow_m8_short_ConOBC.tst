<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m8_short_ConOBC.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m8_short_ConOBC.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.5291e-08</value>
            <value variable="v" tolerance="1e-11">9.82199e-09</value>
            <value variable="p" tolerance="1e-11">1.15756e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">9.53822e-08</value>
            <value variable="v" tolerance="1e-11">5.67855e-08</value>
	    <value variable="p" tolerance="2e-11">2.93157e-07</value>
        </metric>
    </metrics>
</test>
