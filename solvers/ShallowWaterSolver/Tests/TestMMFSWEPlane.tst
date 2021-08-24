<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>MMF SWE solver, DG, P=4</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>TestMMFSWEPlane.xml</parameters>
    <files>
        <file description="Session File">TestMMFSWEPlane.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="eta" tolerance="1e-8">0.000552485</value>
            <value variable="u" tolerance="1e-8">0.000625905</value>
            <value variable="v" tolerance="1e-8">0.0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="eta" tolerance="1e-5">0.031307</value>
            <value variable="u" tolerance="1e-5">0.0233803</value>
            <value variable="v" tolerance="1e-5">0.0</value>
        </metric>
    </metrics>
</test>


