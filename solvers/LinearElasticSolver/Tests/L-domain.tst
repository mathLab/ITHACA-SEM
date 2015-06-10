<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Linear elastic solver P=4</description>
    <executable>LinearElasticSolver</executable>
    <parameters>L-domain.xml</parameters>
    <files>
        <file description="Session File">L-domain.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.000945546</value>
            <value variable="v" tolerance="1e-12">0.00108137</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.014629</value>
            <value variable="v" tolerance="1e-12">0.0234635</value>
        </metric>
    </metrics>
</test>
