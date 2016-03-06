<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Transient Growth (Arpack): Backward-facing step</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>bfs_tg-AR.xml</parameters>
    <files>
        <file description="Session File">bfs_tg-AR.xml</file>
	<file description="Session File">bfs_tg-AR.bse</file>
	<file description="Session File">bfs_tg-AR.rst</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value tolerance="0.001">1.1398,0</value>
        </metric>
    </metrics>
</test>

