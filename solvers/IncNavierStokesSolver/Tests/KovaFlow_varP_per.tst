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
            <value variable="u" tolerance="1e-12">0.000221236</value>
            <value variable="v" tolerance="1e-12">0.000170459</value>
	    <value variable="p" tolerance="1e-12">0.001636</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000366222</value>
            <value variable="v" tolerance="1e-12">0.000395049</value>
	    <value variable="p" tolerance="1e-12">0.0039944</value>
        </metric>
    </metrics>
</test>

