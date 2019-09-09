<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Kovasznay Flow variable P, periodic BC (parallel)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-scotch -I GlobalSysSoln=XxtMultiLevelStaticCond KovaFlow_varP_per.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">KovaFlow_varP_per.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.24657e-05</value>
            <value variable="v" tolerance="1e-12">0.000139872</value>
            <value variable="p" tolerance="1e-12">0.000546763</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000175798</value>
            <value variable="v" tolerance="1e-12">0.000331233</value>
            <value variable="p" tolerance="1e-12">0.00303366</value>
        </metric>
    </metrics>
</test>
