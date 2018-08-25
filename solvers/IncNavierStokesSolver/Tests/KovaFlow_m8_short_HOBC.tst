<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m8_short_HOBC.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m8_short_HOBC.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">1.84877e-08</value>
            <value variable="v" tolerance="1e-11">3.17859e-09</value>
            <value variable="p" tolerance="1e-11">9.68697e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">4.74541e-08</value>
            <value variable="v" tolerance="1e-11">1.17369e-08</value>
	    <value variable="p" tolerance="1e-10">2.92322e-07</value>
        </metric>
    </metrics>
</test>
