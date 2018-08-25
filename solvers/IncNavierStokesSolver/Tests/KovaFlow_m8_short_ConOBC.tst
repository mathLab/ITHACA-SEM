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
            <value variable="u" tolerance="1e-11">1.84566e-08</value>
            <value variable="v" tolerance="1e-11">3.25308e-09</value>
            <value variable="p" tolerance="1e-11">9.98319e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">4.74535e-08</value>
            <value variable="v" tolerance="1e-11">1.20066e-08</value>
	    <value variable="p" tolerance="2e-11">2.92278e-07</value>
        </metric>
    </metrics>
</test>
