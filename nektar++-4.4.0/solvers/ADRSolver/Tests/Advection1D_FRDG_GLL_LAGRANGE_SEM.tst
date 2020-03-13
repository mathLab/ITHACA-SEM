<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady FRDG advection GLL_LAGRANGE_SEM, P=3</description>
    <executable>ADRSolver</executable>
    <parameters>Advection1D_FRDG_GLL_LAGRANGE_SEM.xml</parameters>
    <files>
        <file description="Session File">Advection1D_FRDG_GLL_LAGRANGE_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00988799</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0194621</value>
        </metric>
    </metrics>
</test>
