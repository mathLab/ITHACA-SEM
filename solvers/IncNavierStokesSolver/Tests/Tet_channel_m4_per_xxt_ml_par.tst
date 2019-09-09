<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Tetrahedral elements, P=4, periodic BCs, XxtMultiLevelStaticCond, par(2)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-scotch --verbose -I GlobalSysSoln=XxtMultiLevelStaticCond Tet_channel_m4_per.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Tet_channel_m4_per.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">0</value>
        </metric>
    </metrics>
</test>
