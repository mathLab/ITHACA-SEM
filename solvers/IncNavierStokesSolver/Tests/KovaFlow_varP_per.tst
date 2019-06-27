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
            <value variable="u" tolerance="1e-12">0.00021316</value>
            <value variable="v" tolerance="1e-12">0.000124323</value>
	    <value variable="p" tolerance="1e-12">0.00138562</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00033996</value>
            <value variable="v" tolerance="1e-12">0.000352533</value>
	    <value variable="p" tolerance="1e-12">0.00265904</value>
        </metric>
    </metrics>
</test>

