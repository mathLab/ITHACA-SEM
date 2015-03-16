<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Transient Growth (Modified Arnoldi): Backward-facing step</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>bfs_tg.xml</parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
        <file description="Session File">bfs_tg.bse</file>
        <file description="Session File">bfs_tg.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.0187703</value>
            <value variable="v" tolerance="1e-12">0.019542</value>
            <value variable="p" tolerance="1e-12">0.0292666</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0434972</value>
            <value variable="v" tolerance="1e-12">0.0468062</value>
            <value variable="p" tolerance="1e-12">0.0247981</value>
        </metric>
    </metrics>
</test>
