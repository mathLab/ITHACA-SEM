<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady DG advection, order 1, P=12, periodic bcs</description>
    <executable>ADRSolver</executable>
    <parameters>Advection_m12_DG_periodic.xml</parameters>
    <files>
        <file description="Session File">Advection_m12_DG_periodic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9">1.39014e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">1.32385e-10</value>
        </metric>
    </metrics>
</test>
