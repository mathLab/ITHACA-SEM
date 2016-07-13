<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow variable P, periodic BC (parallel)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-metis -I GlobalSysSoln=XxtMultiLevelStaticCond KovaFlow_varP_per.xml</parameters>
    <processes> 3 </processes>
    <files>
        <file description="Session File">KovaFlow_varP_per.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.74673e-05</value>
            <value variable="v" tolerance="1e-12">0.000201486</value>
	    <value variable="p" tolerance="1e-12">0.000565444</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00011373</value>
            <value variable="v" tolerance="1e-12"> 0.000382783</value>
	    <value variable="p" tolerance="1e-12">0.00344565</value>
        </metric>
    </metrics>
</test>

