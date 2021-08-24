<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady WeakDG advection GAUSS_LAGRANGE, P=3 Q=5</description>
    <executable>ADRSolver</executable>
    <parameters>Advection1D_WeakDG_GAUSS_LAGRANGE.xml</parameters>
    <files>
        <file description="Session File">Advection1D_WeakDG_GAUSS_LAGRANGE.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00961154</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0181693</value>
        </metric>
    </metrics>
</test>
