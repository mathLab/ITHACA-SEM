<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, hexahedra, order 1, P=10,periodic bcs</description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_m10_DG_hex_periodic_nodal.xml</parameters>
    <files>
        <file description="Session File">Advection3D_m10_DG_hex_periodic_nodal.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 3.59823e-05 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-03"> 0.000532131 </value>
        </metric>
    </metrics>
</test>
