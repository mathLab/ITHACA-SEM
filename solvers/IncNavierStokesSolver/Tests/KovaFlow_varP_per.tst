<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow variable P, periodic BC</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-I GlobalSysSoln=DirectMultiLevelStaticCond KovaFlow_varP_per.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_varP_per.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.68743e-05</value>
            <value variable="v" tolerance="1e-12">0.000123899</value>
	    <value variable="p" tolerance="1e-12">0.000569749</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00014278</value>
            <value variable="v" tolerance="1e-12">0.000409628</value>
	    <value variable="p" tolerance="1e-12">0.00342738</value>
        </metric>
    </metrics>
</test>

