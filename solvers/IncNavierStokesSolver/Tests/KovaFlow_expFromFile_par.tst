<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow, expansion from restart file, par(2)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-scotch KovaFlow_expFromFile.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">KovaFlow_expFromFile.xml</file>
	<file description="Session File">KovaFlow_m8.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">3.22767e-07</value>
            <value variable="v" tolerance="1e-12">1.68898e-06</value>
	    <value variable="p" tolerance="1e-12">9.87849e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.31288e-07</value>
            <value variable="v" tolerance="1e-12">2.23918e-06</value>
	    <value variable="p" tolerance="1e-12">3.43342e-05</value>
        </metric>
    </metrics>
</test>
