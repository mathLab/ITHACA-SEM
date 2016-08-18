<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>ADRSolver C0 projection</description>
    <executable>ADRSolver</executable>
    <parameters>Projection2D.xml</parameters>
    <files>
        <file description="Session File">Projection2D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="0">
            <value variable="u" tolerance="1e-12">1.19001e-16</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="1e-12">4.16334e-16</value>
        </metric>
    </metrics>
</test>
