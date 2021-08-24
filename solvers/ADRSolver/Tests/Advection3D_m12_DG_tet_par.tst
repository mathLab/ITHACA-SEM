<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, tetrahedra, order 4, P=12, par(3)</description>
    <executable>ADRSolver</executable>
    <parameters>--use-scotch Advection3D_m12_DG_tet.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Advection3D_m12_DG_tet.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.00034e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">5.92099e-06</value>
        </metric>
    </metrics>
</test>



