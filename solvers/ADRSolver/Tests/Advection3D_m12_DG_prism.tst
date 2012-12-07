<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, prisms, order 4, P=14</description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_m12_DG_prism.xml</parameters>
    <files>
        <file description="Session File">Advection3D_m12_DG_prism.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">6.72764e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.01404e-05</value>
        </metric>
    </metrics>
</test>


