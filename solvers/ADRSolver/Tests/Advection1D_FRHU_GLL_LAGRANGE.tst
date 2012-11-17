<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady FRHU advection GLL_LAGRANGE, P=3</description>
    <executable>ADRSolver</executable>
    <parameters>Advection1D_FRHU_GLL_LAGRANGE.xml</parameters>
    <files>
        <file description="Session File">Advection1D_FRHU_GLL_LAGRANGE.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.0369431</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0642327</value>
        </metric>
    </metrics>
</test>
