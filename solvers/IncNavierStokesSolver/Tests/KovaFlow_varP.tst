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
            <value variable="u" tolerance="1e-12">9.42175e-05</value>
            <value variable="v" tolerance="1e-12">7.41865e-06</value>
	    <value variable="p" tolerance="2e-9">0.000981868</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000152464</value>
            <value variable="v" tolerance="1e-12">9.02121e-06</value>
	    <value variable="p" tolerance="1e-12">0.00296616</value>
        </metric>
    </metrics>
</test>
