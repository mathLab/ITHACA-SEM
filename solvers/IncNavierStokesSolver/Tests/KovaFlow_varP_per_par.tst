<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Kovasznay Flow variable P, periodic BC (parallel)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-metis -I GlobalSysSoln=XxtMultiLevelStaticCond KovaFlow_varP_per.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">KovaFlow_varP_per.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.18634e-05</value>
            <value variable="v" tolerance="1e-12">0.000139364</value>
            <value variable="p" tolerance="1e-12">0.000546725</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000174552</value>
            <value variable="v" tolerance="1e-12">0.000331172</value>
            <value variable="p" tolerance="1e-12">0.00303379</value>
        </metric>
    </metrics>
</test>
