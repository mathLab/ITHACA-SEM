<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, tetrahedra, order 4, P=Variable</description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_m10_DG_hex_VarP.xml</parameters>
    <files>
        <file description="Session File">Advection3D_m10_DG_hex_VarP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.67026e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">1.94197e-05</value>
        </metric>
    </metrics>
</test>
