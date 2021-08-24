<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=8 Weak Pressure VSS</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m8_short_HOBC_VCSWeakPress.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m8_short_HOBC_VCSWeakPress.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.51879e-08</value>
            <value variable="v" tolerance="1e-11">9.55783e-09</value>
	    <value variable="p" tolerance="1e-11">1.11241e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">9.42611e-08</value>
            <value variable="v" tolerance="1e-11">5.59167e-08</value>
	    <value variable="p" tolerance="1e-10">2.94001e-07</value>
        </metric>
    </metrics>
</test>
