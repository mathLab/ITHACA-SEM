<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady FRDG advection GLL_LAGRANGE, P=3</description>
    <executable>ADRSolver</executable>
    <parameters>Advection1D_FRDG_GLL_LAGRANGE.xml</parameters>
    <files>
        <file description="Session File">Advection1D_FRDG_GLL_LAGRANGE.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00960004</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0177832</value>
        </metric>
    </metrics>
</test>
