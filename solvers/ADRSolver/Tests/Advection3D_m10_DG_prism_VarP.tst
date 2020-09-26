<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, tetrahedra, order 4, P=Variable</description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_m10_DG_prism_VarP.xml</parameters>
    <files>
        <file description="Session File">Advection3D_m10_DG_prism_VarP.xml</file>
    </files>
    <metrics>
            <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.5189e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">7.42306e-06</value>
        </metric>
    </metrics>
</test>

