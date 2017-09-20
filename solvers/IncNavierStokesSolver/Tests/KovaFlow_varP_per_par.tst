<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Kovasznay Flow variable P, periodic BC (parallel)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-metis -I GlobalSysSoln=XxtMultiLevelStaticCond KovaFlow_varP_per.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">KovaFlow_varP_per.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">7.50905e-05</value>
            <value variable="v" tolerance="1e-12">0.000120802</value>
            <value variable="p" tolerance="1e-12">0.000545616</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000142621</value>
            <value variable="v" tolerance="1e-12">0.000329214</value>
            <value variable="p" tolerance="1e-12">0.00303305</value>
        </metric>
    </metrics>
</test>
