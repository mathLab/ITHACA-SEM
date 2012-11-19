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
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.0091919</value>
            <value variable="v" tolerance="1e-12">0.00706682</value>
            <value variable="p" tolerance="1e-12">0.191912</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">0</value>
        </metric>
    </metrics>
</test>

