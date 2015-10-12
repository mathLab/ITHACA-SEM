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
        <metric type="Eigenvalue" id="0">
            <value tolerance="0.001">1.1398,0</value>
        </metric>
    </metrics>
</test>
