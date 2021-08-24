<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady FRSD inviscidBurger GLL_LAGRANGE_SEM, P=10</description>
    <executable>ADRSolver</executable>
    <parameters>InviscidBurger1D_FRSD_GLL_LAGRANGE_SEM.xml</parameters>
    <files>
        <file description="Session File">InviscidBurger1D_FRSD_GLL_LAGRANGE_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">2.28218</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.5</value>
        </metric>
    </metrics>
</test>
