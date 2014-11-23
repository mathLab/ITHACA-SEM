<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Transient Growth (Modified Arnoldi): Backward-facing step</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>bfs_tg.xml</parameters>
    <processes> 2 </processes>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.bse</file>
	<file description="Session File">bfs_tg.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.0187569</value>
            <value variable="v" tolerance="1e-12">0.0195275</value>
            <value variable="p" tolerance="1e-6">0.0171521</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0434466</value>
            <value variable="v" tolerance="1e-12">0.0467714</value>
            <value variable="p" tolerance="1e-6">0.0219808</value>
        </metric>
    </metrics>
</test>

